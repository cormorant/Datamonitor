#!/usr/bin/env python
# coding: utf-8
from __future__ import division
"""
Описание формата Байкал-5
"""
APP_NAME = "BaikalLib"
__version__="0.0.1.2"
COMPANY_NAME = 'GIN SB RAS'

import os
import struct
import ctypes
import datetime

import numpy as np
#from . import sigtools
import sigtools

#===============================================================================
class GeneralHeader53(ctypes.LittleEndianStructure):
    _fields_ = [("number_of_channel", ctypes.c_uint16),
                ("test_type", ctypes.c_uint16),
                ("version", ctypes.c_uint16),
                ("day", ctypes.c_uint16),
                ("month", ctypes.c_uint16),
                ("year", ctypes.c_uint16),
                ("satellites", ctypes.c_uint16),
                ("end_invalid", ctypes.c_uint16),
                ("synchronize", ctypes.c_uint16),
                ("digits", ctypes.c_uint16),
                ("start_invalid", ctypes.c_uint16),
                ("constant", ctypes.c_uint16),
                ("acq_version", ctypes.c_uint16),
                ("system_freq", ctypes.c_uint16),
                ("max_criteria", ctypes.c_uint16),
                ("satellite_instart", ctypes.c_uint16),
                ("station_name", ctypes.c_char * 16),
                ("sample_period", ctypes.c_double),
                ("time_begin", ctypes.c_double),
                ("correction", ctypes.c_double),
                ("latitude", ctypes.c_double),
                ("longitude", ctypes.c_double),
                ("before_synch", ctypes.c_uint64),
                ("synch_point", ctypes.c_uint64),
                ("after_synch", ctypes.c_uint64),
                ("synch_point_start", ctypes.c_uint32),
                ("reserved", ctypes.c_uint32)]

class ChannelHeader(ctypes.LittleEndianStructure):
    _fields_ = [("channel_number", ctypes.c_short),
                ("reserved", ctypes.c_char * 6),
                ("channel_name", ctypes.c_char * 24),
                ("channel_mark", ctypes.c_char * 24),
                ("channel_rate", ctypes.c_double),
                ("reserved2", ctypes.c_double)]


def _is_xx(filename):
    """ Checks whether a file is Baykal XX waveform data or not """
    with open(filename, "rb") as _f:
        try:
            number_of_channels, = struct.unpack("h", _f.read(2))
        except struct.error:
            return False
    # должно быть вразумительное число каналов
    if 0 < number_of_channels <= 12: return True


def _read_xx53_channel(filename, filter=False):
    """ Reads an Baykal XX v53 waveform file and returns a Stream object """
    general = GeneralHeader53()
    headers = []
    traces = []
    with open(filename, 'rb') as _f:
        _f.readinto(general)
        number_of_channels = general.number_of_channel
        # check year
        year = general.year
        if year < 1900: year += 2000
        # get start time of traces
        starttime = datetime.datetime(year, general.month, general.day)
        header = {}
        header['sampling_rate'] = 1. / general.sample_period #TODO: add round...
        # use only 3 first symbols
        header['station'] = general.station_name[:3].upper()# may have 4 symbols
        # t0 - secodns from 00:00:00 of day
        delta = datetime.timedelta(seconds=general.time_begin)
        header['starttime'] = starttime + delta
        # read all channel headers (for what?)
        #for _num in range(number_of_channels):
        #    ch = ChannelHeader()
        #    _f.readinto(ch)
        #    headers.append(ch)
        # read data section
        data = np.frombuffer(_f.read(-1),
            dtype=np.int16 if general.digits==16 else np.int32)
        # get data as float
        data = data.astype(np.float32)
    #=== take all channels, demultiplexed
    # обрезать массив с конца пока он не делится на 3
    while data.size % 3 != 0: data = data[:-1]
    # demultiplexing
    len_of_channel = int(data.size / number_of_channels)
    # reshape data
    data = data.reshape(len_of_channel, number_of_channels)
    data = data.T
    # return data and header
    return data, header

#===============================================================================
# add filtering functionality here...

#from scipy.signal import butter, lfilter#, freqz
#================================================

def zpk2tf(z, p, k):
    """ Return polynomial transfer function representation from zeros and poles """
    z = np.atleast_1d(z)
    k = np.atleast_1d(k)
    if len(z.shape) > 1:
        temp = np.poly(z[0])
        b = np.zeros((z.shape[0], z.shape[1] + 1), temp.dtype.char)
        if len(k) == 1:
            k = [k[0]] * z.shape[0]
        for i in range(z.shape[0]):
            b[i] = k[i] * np.poly(z[i])
    else:
        b = k * np.poly(z)
    a = np.atleast_1d(np.poly(p))
    # Use real output if possible
    if issubclass(b.dtype.type, np.complexfloating):
        # if complex roots are all complex conjugates, the roots are real.
        roots = np.asarray(z, complex)
        pos_roots = np.compress(roots.imag > 0, roots)
        neg_roots = np.conjugate(np.compress(roots.imag < 0, roots))
        if len(pos_roots) == len(neg_roots):
            if np.all(np.sort_complex(neg_roots) ==
                         np.sort_complex(pos_roots)):
                b = b.real.copy()

    if issubclass(a.dtype.type, np.complexfloating):
        # if complex roots are all complex conjugates, the roots are real.
        roots = np.asarray(p, complex)
        pos_roots = np.compress(roots.imag > 0, roots)
        neg_roots = np.conjugate(np.compress(roots.imag < 0, roots))
        if len(pos_roots) == len(neg_roots):
            if np.all(np.sort_complex(neg_roots) ==
                         np.sort_complex(pos_roots)):
                a = a.real.copy()
    return b, a

def _relative_degree(z, p):
    """ Return relative degree of transfer function from zeros and poles """
    degree = len(p) - len(z)
    if degree < 0: raise ValueError("Improper transfer function.(((")
    else: return degree

def _zpklp2lp(z, p, k, wo=1.0):
    """ Transform a lowpass filter prototype to a different frequency """
    z = np.atleast_1d(z)
    p = np.atleast_1d(p)
    wo = float(wo)  # Avoid int wraparound
    degree = _relative_degree(z, p)
    # Scale all points radially from origin to shift cutoff frequency
    z_lp = wo * z
    p_lp = wo * p
    # Each shifted pole decreases gain by wo, each shifted zero increases it.
    # Cancel out the net change to keep overall gain the same
    k_lp = k * wo**degree
    return z_lp, p_lp, k_lp

def _zpkbilinear(z, p, k, fs):
    """ Return a digital filter from an analog one using a bilinear transform """
    z = np.atleast_1d(z)
    p = np.atleast_1d(p)
    degree = _relative_degree(z, p)
    fs2 = 2.0 * fs
    # Bilinear transform the poles and zeros
    z_z = (fs2 + z) / (fs2 - z)
    p_z = (fs2 + p) / (fs2 - p)
    # Any zeros that were at infinity get moved to the Nyquist frequency
    z_z = np.append(z_z, -np.ones(degree))
    # Compensate for gain change
    k_z = k * np.real(np.prod(fs2 - z) / np.prod(fs2 - p))
    return z_z, p_z, k_z

def buttap(N):
    """ Return (z,p,k) for analog prototype of Nth-order Butterworth filter """
    if abs(int(N)) != N:
        raise ValueError("Filter order must be a nonnegative integer")
    z = np.array([])
    m = np.arange(-N+1, N, 2)
    # Middle value is 0 to ensure an exactly real pole
    p = -np.exp(1j * np.pi * m / (2 * N))
    k = 1
    return z, p, k

'''
def buttord(wp, ws, gpass, gstop, analog=False):
    """ Butterworth filter order selection """
    wp = np.atleast_1d(wp)
    ws = np.atleast_1d(ws)
    filter_type = 2 * (len(wp) - 1)
    filter_type += 1
    if wp[0] >= ws[0]:
        filter_type += 1
    # Pre-warp frequencies for digital filter design
    if not analog:
        passb = np.tan(np.pi * wp / 2.0)
        stopb = np.tan(np.pi * ws / 2.0)
    else:
        passb = wp * 1.0
        stopb = ws * 1.0

    if filter_type == 1:            # low
        nat = stopb / passb
    elif filter_type == 2:          # high
        nat = passb / stopb
    elif filter_type == 3:          # stop
        wp0 = optimize.fminbound(band_stop_obj, passb[0], stopb[0] - 1e-12,
                                 args=(0, passb, stopb, gpass, gstop,
                                       'butter'),
                                 disp=0)
        passb[0] = wp0
        wp1 = optimize.fminbound(band_stop_obj, stopb[1] + 1e-12, passb[1],
                                 args=(1, passb, stopb, gpass, gstop,
                                       'butter'),
                                 disp=0)
        passb[1] = wp1
        nat = ((stopb * (passb[0] - passb[1])) /
               (stopb ** 2 - passb[0] * passb[1]))
    elif filter_type == 4:          # pass
        nat = ((stopb ** 2 - passb[0] * passb[1]) /
               (stopb * (passb[0] - passb[1])))

    nat = min(abs(nat))

    GSTOP = 10 ** (0.1 * abs(gstop))
    GPASS = 10 ** (0.1 * abs(gpass))
    ord = int(ceil(log10((GSTOP - 1.0) / (GPASS - 1.0)) / (2 * log10(nat))))

    # Find the Butterworth natural frequency WN (or the "3dB" frequency")
    # to give exactly gpass at passb.
    try:
        W0 = (GPASS - 1.0) ** (-1.0 / (2.0 * ord))
    except ZeroDivisionError:
        W0 = 1.0
        print("Warning, order is zero...check input parameters.")

    # now convert this frequency back from lowpass prototype
    # to the original analog filter

    if filter_type == 1:  # low
        WN = W0 * passb
    elif filter_type == 2:  # high
        WN = passb / W0
    elif filter_type == 3:  # stop
        WN = numpy.zeros(2, float)
        discr = sqrt((passb[1] - passb[0]) ** 2 +
                     4 * W0 ** 2 * passb[0] * passb[1])
        WN[0] = ((passb[1] - passb[0]) + discr) / (2 * W0)
        WN[1] = ((passb[1] - passb[0]) - discr) / (2 * W0)
        WN = numpy.sort(abs(WN))
    elif filter_type == 4:  # pass
        W0 = numpy.array([-W0, W0], float)
        WN = (-W0 * (passb[1] - passb[0]) / 2.0 +
              sqrt(W0 ** 2 / 4.0 * (passb[1] - passb[0]) ** 2 +
                   passb[0] * passb[1]))
        WN = numpy.sort(abs(WN))
    else:
        raise ValueError("Bad type: %s" % filter_type)

    if not analog:
        wn = (2.0 / pi) * arctan(WN)
    else:
        wn = WN

    if len(wn) == 1:
        wn = wn[0]
    return ord, wn
'''

def iirfilter(N, Wn, rp=None, rs=None, btype='band', analog=False,
              ftype='butter', output='ba'):
    """ IIR digital and analog filter design given order and critical points """
    ftype, btype, output = [x.lower() for x in (ftype, btype, output)]
    Wn = np.asarray(Wn)
    try:
        btype = band_dict[btype]
    except KeyError:
        raise ValueError("'%s' is an invalid bandtype for filter." % btype)

    try:
        typefunc = filter_dict[ftype][0]
    except KeyError:
        raise ValueError("'%s' is not a valid basic IIR filter." % ftype)

    if output not in ['ba', 'zpk', 'sos']:
        raise ValueError("'%s' is not a valid output form." % output)

    if rp is not None and rp < 0:
        raise ValueError("passband ripple (rp) must be positive")

    if rs is not None and rs < 0:
        raise ValueError("stopband attenuation (rs) must be positive")
    # Get analog lowpass prototype
    if typefunc == buttap:
        z, p, k = typefunc(N)
    else:
        raise NotImplementedError("'%s' not implemented in iirfilter." % ftype)
    # Pre-warp frequencies for digital filter design
    if not analog:
        if np.any(Wn <= 0) or np.any(Wn >= 1):
            raise ValueError("Digital filter critical freqs must be 0< Wn <1")
        fs = 2.0
        warped = 2 * fs * np.tan(np.pi * Wn / fs)
    else:
        warped = Wn
    # transform to lowpass, bandpass, highpass, or bandstop
    if btype in ('lowpass', 'highpass'):
        if np.size(Wn) != 1:
            raise ValueError('Must specify a single critical frequency Wn')
        if btype == 'lowpass':
            z, p, k = _zpklp2lp(z, p, k, wo=warped)
        #elif btype == 'highpass':
        #    z, p, k = _zpklp2hp(z, p, k, wo=warped)
    else:
        raise NotImplementedError("'%s' not implemented in iirfilter." % btype)
    # Find discrete equivalent if necessary
    if not analog:
        z, p, k = _zpkbilinear(z, p, k, fs=fs)

    # Transform to proper out type (pole-zero, state-space, numer-denom)
    if output == 'zpk':
        return z, p, k
    elif output == 'ba':
        return zpk2tf(z, p, k)
    #elif output == 'sos':
    #    return zpk2sos(z, p, k)


def butter(N, Wn, btype='low', analog=False, output='ba'):
    """ Butterworth digital and analog filter design """
    return iirfilter(N, Wn, btype=btype, analog=analog,
                     output=output, ftype='butter')

#=== lfilter
def lfilter(b, a, x, axis=-1, zi=None):
    """ Filter data along one-dimension with an IIR or FIR filter """
    a = np.atleast_1d(a)
    if len(a) == 1:
        b = np.asarray(b)
        a = np.asarray(a)
        if b.ndim != 1 and a.ndim != 1:
            raise ValueError('object of too small depth for desired array')
        x = np.asarray(x)
        inputs = [b, a, x]
        if zi is not None:
            # _linear_filter does not broadcast zi, but does do expansion of singleton dims
            zi = np.asarray(zi)
            if zi.ndim != x.ndim:
                raise ValueError('object of too small depth for desired array')
            expected_shape = list(x.shape)
            expected_shape[axis] = b.shape[0] - 1
            expected_shape = tuple(expected_shape)
            # check the trivial case where zi is the right shape first
            if zi.shape != expected_shape:
                strides = zi.ndim * [None]
                if axis < 0:
                    axis += zi.ndim
                for k in range(zi.ndim):
                    if k == axis and zi.shape[k] == expected_shape[k]:
                        strides[k] = zi.strides[k]
                    elif k != axis and zi.shape[k] == expected_shape[k]:
                        strides[k] = zi.strides[k]
                    elif k != axis and zi.shape[k] == 1:
                        strides[k] = 0
                    else:
                        raise ValueError('Unexpected shape for zi: expected '
                                         '%s, found %s.' %
                                         (expected_shape, zi.shape))
                zi = np.lib.stride_tricks.as_strided(zi, expected_shape, strides)
            inputs.append(zi)
        dtype = np.result_type(*inputs)

        if dtype.char not in 'fdgFDGO':
            raise NotImplementedError("input type '%s' not supported" % dtype)

        b = np.array(b, dtype=dtype)
        a = np.array(a, dtype=dtype, copy=False)
        b /= a[0]
        x = np.array(x, dtype=dtype, copy=False)

        out_full = np.apply_along_axis(lambda y: np.convolve(b, y), axis, x)
        ind = out_full.ndim * [slice(None)]
        if zi is not None:
            ind[axis] = slice(zi.shape[axis])
            out_full[ind] += zi

        ind[axis] = slice(out_full.shape[axis] - len(b) + 1)
        out = out_full[ind]

        if zi is None:
            return out
        else:
            ind[axis] = slice(out_full.shape[axis] - len(b) + 1, None)
            zf = out_full[ind]
            return out, zf
    else:
        if zi is None:
            return sigtools._linear_filter(b, a, x, axis)
        else:
            return sigtools._linear_filter(b, a, x, axis, zi)
#=== lfilter ends


filter_dict = {
    'butter': [buttap, 'buttord'],
    'butterworth': [buttap, 'buttord'],
}

band_dict = {
    'band': 'bandpass',
    'bandpass': 'bandpass',
    'pass': 'bandpass',
    'bp': 'bandpass',

    'bs': 'bandstop',
    'bandstop': 'bandstop',
    'bands': 'bandstop',
    'stop': 'bandstop',

    'l': 'lowpass',
    'low': 'lowpass',
    'lowpass': 'lowpass',
    'lp': 'lowpass',

    #'high': 'highpass',
    #'highpass': 'highpass',
    #'h': 'highpass',
    #'hp': 'highpass',
}
#================================================


def butter_lowpass(cutoff, fs, order):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, fs, cutoff=10, order=4):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y
