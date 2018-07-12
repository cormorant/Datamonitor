#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  windatamonitor.py
#  
"""
    Мониторинг изменений в папке, куда сохраняется поток данных с регистратора
Ангара-7Б или Байкал-11, или (в перспективе) Иркут-24 (через фтп, wget).

    Algorithm STA/LTA and marking file works correctly.
But now lets make separate program that uses couple (1, 2 or 3)
of  XX files as input and writes resulting 3-minutes (usually) file'

    The default filter is RMHP(10)»ITAPER(30)»BW(3,0.7,2.0)»STALTA(2,80).
Recommendations for other distance ranges are:
I teleseismic: BW(3,0.7,2.0)»STALTA(2,80)
I regional: BW(4,1.0,10.0)»STALTA(1,40)
I local: BW(4,4.0,20.0)»STALTA(0.2,10)

    Lets also add filter, cause we can have technological high-freq signals.

    Писать информацию в лог-файл, когда появился ФАЙЛ и информацию:
    DATETIME, FILENAME, START-OF-FILE, END-OF-FILE, + info (MAX, votes)

    TODO: add verbose mode.
"""
from __future__ import division
import os
import sys
import datetime; now = datetime.datetime.now
import time
import shutil

from watchdog.observers import Observer
from watchdog.events import (FileSystemEventHandler,
    DirModifiedEvent, FileModifiedEvent)

import ConfigParser

from baikal import _is_xx, _read_xx53_channel, butter_lowpass_filter

from trigger import classic_sta_lta, trigger_onset


# папка откуда запускается программа
def module_path():
    if hasattr(sys, "frozen"):
        return os.path.dirname(sys.executable)
    return os.path.dirname(__file__)
# get Current Dir
CurrDir = module_path()

# if not LOG dir - create it!
LOG_DIR = os.path.join(CurrDir, "log")
if not os.path.exists(LOG_DIR):
    os.mkdir(LOG_DIR)
sys.stderr = open(os.path.join(CurrDir, "log", "windatamonitor.err"), "a")

# register file - write there datetime and filename that was triggered
#REGISTER_FILE = open(os.path.join(CurrDir, "log", "triggered.register"), "a")


DEFAULT_SETTINGS = {
    'INTERVAL': 20,
    'MAX_LEVEL': 10.,
    'MINIMAL_LENGTH': 4,
    'MIN_LEVEL': 0.5,
    'SAMPLES': 6000,
    'PATH': 'c:/angara/files', # default path for angara files
    'WINDOW_LONG': 30,
    'WINDOW_SHORT': 1,
    "CHANNEL": 2,
    # filter or not
    "FILTER": False,
    # watch just 1 channel
    "WATCHING": 1,
    # seconds
    "SECONDS_BEFORE": 90,
    "SECONDS_AFTER": 90,
}


def read_configfile(configfilename, section="main"):
    """ config file (settings.conf) reading; общие настройки """
    Settings = {}
    # считывание настроек
    config = ConfigParser.SafeConfigParser()
    config.read(configfilename)
    # path
    Settings["PATH"] = config.get(section, 'path')
    # read INT values
    KEYS = ("SAMPLES", "INTERVAL", "WINDOW_SHORT", "WINDOW_LONG",
        "MINIMAL_LENGTH", "SECONDS_BEFORE", "SECONDS_AFTER", "WATCHING")
    for key in KEYS:
        Settings[key] = config.getint(section, key)
    # get FLOAT values
    for key in ("MAX_LEVEL", "MIN_LEVEL", "CUTOFF"):
        Settings[key] = config.getfloat(section, key)
    # boolean value
    key = "FILTER"
    Settings[key] = config.getboolean(section, key)
    return Settings


class MyHandler(FileSystemEventHandler):
    """ FileSystemEventHandler - базовый класс обработчика событий изменения """
    def __init__(self):
        super(MyHandler, self).__init__()

    def on_modified(self, event):
        """ runs when some file in Dir is created or midified """
        # here may be DirModifiedEvent, we look for FileModifiedEvent
        if not isinstance(event, DirModifiedEvent):
            #===================================================================
            try:
                #=== real work here: what happends when we have new file!!!!!!!!
                self.trigger(event.src_path)
                #===============================================================
            except BaseException as e:
                print(e)
        #=======================================================================


    def my_rolling_window(self, arrays, window, yield_last_array=True):
        """ Native rolling_window may cause problem with memory.
        Thats it, lets use just an ordinal cycle For... """
        # calc step: usually 3 minutes, so window / 3 = 1 minute
        step = int(window / 3)
        # length of all array
        length_array = arrays[0].size
        # start and end indices
        start_index = 0
        while True:
            # calc end of array
            end_index = start_index + window
            if end_index > length_array: break
            yield arrays[:, start_index:end_index]
            # increase start_index with step
            start_index += step
        # finally, yield last array, from the end
        if yield_last_array:
            yield arrays[:, -window:]

    def trigger(self, filename):
        """ analyze file minute by minute """
        # check if it XX file
        if not _is_xx(filename): return
        #===
        print("\n+++\t%s" % filename)
        # read file
        result = _read_xx53_channel(filename)
        if result is None:
            print("\nProblem occured while reading file %s!" % filename)
            return
        else:
            # if valid file, there must be data array and header dict
            data, header = result
        #===
        #= if trace is OK...
        df = header['sampling_rate'] # sample rate (in Hz)
        # first of all, check length of returned array: must be > N
        TOTAL_SAMPLES = Settings["SAMPLES"]
        # on first channel, check length...
        data_length = data[0].size
        if data_length < TOTAL_SAMPLES:
            print("\nNot enough samples in file! Lets wait...")
            return
        #=== analyzing starts... But have to use all 3 channels!
        starttime = header['starttime'] # start time of trace
        # calc end time of data in this file
        #endtime = starttime + datetime.timedelta(seconds=data_length / df)
        #===
        # short and long window
        nsta, nlta = int(WINDOW_SHORT * df), int(WINDOW_LONG * df)
        # скользящее (?) окно для анализа файла (по 3 минуты обычно)
        for ns, ew, z in self.my_rolling_window(data, TOTAL_SAMPLES):
            #=== lets run trigger to find Event
            # calculating the characteristic function - для всех каналов
            sum_votes = 0
            valid_on_off = None
            for a in (ns, ew, z):
                # before calculating, lets use filtering, if needed
                if Settings["FILTER"]:
                    a = butter_lowpass_filter(a, header["sampling_rate"],
                        cutoff=Settings["CUTOFF"])
                # calc characterizing function (STA/LTA)
                cft = classic_sta_lta(a, nsta, nlta)
                print cft.max()
                # получить номер отсчета с которого включился триггер
                on_off = self.processCFT(cft)
                if on_off is not None:
                    # сначала проверить сколько секунд между ВКЛ и ВЫКЛ
                    # номера отсчетов где триггер ВКЛ и ВЫКЛ
                    SamplesOn, SamplesOff = on_off[0]
                    # считать сколько это в секундах
                    gap_in_sec = (SamplesOff - SamplesOn) / df
                    print "GAP:", gap_in_sec
                    # промежуток в СЕК должен быть больше МИНИМУМА (4 сек)
                    # но если амплитуды о-очень большие (*2 >), достаточно 1 сек
                    if (gap_in_sec >= MINIMAL_LENGTH) or ((gap_in_sec >= 1) and
                        (cft.max() >= MAX_LEVEL * 2)):
                        sum_votes += 1
                        valid_on_off = on_off
            #===========================================================
            # сработать должны для 2х каналов из 3х (см. настройки)
            print "Votes:", sum_votes
            if sum_votes >= Settings["WATCHING"]:
                #================
                print("!"*33)
                print("\tRun trigger on event (max=%.1f, votes=%d)." %
                    (cft.max(), sum_votes))
                self.trigger_is_on(filename, valid_on_off)
                print("!"*33)
                #================
                # finally break moving window...
                # we don't want multiple triggers in one file
                break
        #print("\n\t...Done!")
        print

    
    def trigger_is_on(self, filename, on_off):
        """ триггер включился... """
        # может просто пометить файл знаком '!' ?
        # остальную работу предоставить SendAgent'у?
        # mark found file with !
        if not "!" in filename:# may be there already ! in file name
            print("Mark file %s with '!'!!!" % filename)
            shutil.move(filename, filename+"!")
        #===
        '''
        # вычислим время включения триггера по количеству отсчетов
        seconds_on = int( round(SamplesOn / df) )
        seconds_off = int( round(SamplesOff / df) )
        TimeTriggerOn = starttime + datetime.timedelta(seconds=seconds_on)
        TimeTriggerOff = starttime + datetime.timedelta(seconds=seconds_off)
        #= increase END and decrease START time by values from config
        time1 = TimeTriggerOn - datetime.timedelta(seconds=SECONDS_BEFORE)
        time2 = TimeTriggerOff + datetime.timedelta(seconds=SECONDS_AFTER)
        print("\n")
        print("!"*33)
        print('Cut from %s to %s' % (time1, time2))
        print("!"*33)
        # проверить, времена вступлений в нашем файле? or wait next
        if time1 < starttime:
            print("\nHave to find PREVIOUS file, to get start...")
        if time2 > endtime:
            print("\nHave to wait for NEXT file, to get the end...")
            msg = "Have to wait for NEXT file, but not Implemented!"
            raise NotImplementedError, msg
        #============================================
        # output info in Register file
        #s = "{0}\t{1}; MAX={2:.2f}; Sec={3:.2f}".format(now(),
        #    filename, cft.max(), gap_in_sec)
        #REGISTER_FILE.write("\n"+s)
        '''

    def processCFT(self, cft):
        """ обработка характеризующей функции """
        # включился ли триггер?
        MAX_CFT = cft.max()
        if MAX_CFT > MAX_LEVEL:
            # calc on and off time
            on_off = trigger_onset(cft, MAX_LEVEL, MIN_LEVEL)
            return on_off
    
def main(Settings):
    """ main func """
    # run observing on FS
    observer = Observer()
    observer.schedule(MyHandler(), path=Settings["PATH"], recursive=True)
    # после вызова start() мы получаем фоновый поток, следящий за изменениями ФС
    observer.start()
    # костыль
    number_of_runs = 0
    try:
        while True:
            # nice output
            s = "\t#%d\tPress <Ctrl+C> to exit..." % number_of_runs
            sys.stdout.write("\r" + s)
            sys.stdout.flush()
            # stop plot
            time.sleep(Settings['INTERVAL'])
            number_of_runs += 1
    # Ждем событий изменений ФС до прихода Ctrl+C (SIGINT)
    except KeyboardInterrupt:
        print("\n")
        print("*"*77)
        print("\tExit...")
        observer.stop()
    observer.join()
    return 0


if __name__ == '__main__':
    # аргумент (файл с настройками)
    if len(sys.argv) > 1:
        configfilename = sys.argv[1]
        if not os.path.exists(configfilename):
            print("Config file '%s' not found! Use Default settings..." %
                configfilename)
        try:
            # read settings
            Settings = read_configfile(configfilename)
            # если нет такого файла - ахтунг!!!
        except BaseException as e:
            print("Error reading config file: %s" % e)
            #sys.exit(1)
            Settings = DEFAULT_SETTINGS
    else:
        # если не указан - использовать настройки по умолчанию
        Settings = DEFAULT_SETTINGS
    # get params from Settings
    MAX_LEVEL, MIN_LEVEL = Settings['MAX_LEVEL'], Settings['MIN_LEVEL']
    # windows...
    WINDOW_SHORT, WINDOW_LONG, MINIMAL_LENGTH = [Settings[k] for k in
        ("WINDOW_SHORT", "WINDOW_LONG", "MINIMAL_LENGTH")]
    # seconds...
    SECONDS_BEFORE, SECONDS_AFTER = (Settings["SECONDS_BEFORE"],
        Settings["SECONDS_AFTER"])
    #===========================================================================
    # also check PATH settings
    path = Settings["PATH"]
    # make full path if not
    if not os.path.isabs(path): path = os.path.join(CurrDir, path)
    if not os.path.exists(path):
        print("Path %s not found! Exit..." % path)
        sys.exit(1)
    #=== Запуск программы
    print("*"*77)
    print("\tProgram started at %s" % now())
    print("*"*77)
    sys.exit(main(Settings))
    #===========================================================================
