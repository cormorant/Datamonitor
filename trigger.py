# -*- coding: utf-8 -*-
from __future__ import division

import ctypes
import numpy as np
from collections import deque


head_stalta_t = np.dtype([
    (str('N'), np.uint32, 1),
    (str('nsta'), np.uint32, 1),
    (str('nlta'), np.uint32, 1),
], align=True)


# Import shared libsignal
_library = 'libsignal.pyd'
cdll = ctypes.CDLL(_library)


cdll.stalta.argtypes = [
    np.ctypeslib.ndpointer(dtype=head_stalta_t, ndim=1,
                           flags=str('C_CONTIGUOUS')),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,
                           flags=str('C_CONTIGUOUS')),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,
                           flags=str('C_CONTIGUOUS')),
]
cdll.stalta.restype = ctypes.c_int



def classic_sta_lta(a, nsta, nlta):
    """ Computes the standard STA/LTA from a given input array a. The length of
    the STA is given by nsta in samples, respectively is the length of the
    LTA given by nlta in samples """
    data = a
    # initialize C struct / NumPy structured array
    head = np.empty(1, dtype=head_stalta_t)
    head[:] = (len(data), nsta, nlta)
    # ensure correct type and contiguous of data
    data = np.ascontiguousarray(data, dtype=np.float64)
    # all memory should be allocated by python
    charfct = np.empty(len(data), dtype=np.float64)
    # run and check the error-code
    errcode = cdll.stalta(head, data, charfct)
    if errcode != 0:
        raise Exception('ERROR %d stalta: len(data) < nlta' % errcode)
    return charfct


def trigger_onset(charfct, thres1, thres2, max_len=9e99, max_len_delete=False):
    """ Calculate trigger on and off times """
    # 1) find indices of samples greater than threshold
    # 2) calculate trigger "of" times by the gap in trigger indices
    #    above the threshold i.e. the difference of two following indices
    #    in ind is greater than 1
    # 3) in principle the same as for "of" just add one to the index to get
    #    start times, this operation is not supported on the compact
    #    syntax
    # 4) as long as there is a on time greater than the actual of time find
    #    trigger on states which are greater than last of state an the
    #    corresponding of state which is greater than current on state
    # 5) if the signal stays above thres2 longer than max_len an event
    #    is triggered and following a new event can be triggered as soon as
    #    the signal is above thres1
    ind1 = np.where(charfct > thres1)[0]
    if len(ind1) == 0:
        return []
    ind2 = np.where(charfct > thres2)[0]
    #
    on = deque([ind1[0]])
    of = deque([-1])
    # determine the indices where charfct falls below off-threshold
    ind2_ = np.empty_like(ind2, dtype=bool)
    ind2_[:-1] = np.diff(ind2) > 1
    # last occurence is missed by the diff, add it manually
    ind2_[-1] = True
    of.extend(ind2[ind2_].tolist())
    on.extend(ind1[np.where(np.diff(ind1) > 1)[0] + 1].tolist())
    # include last pick if trigger is on or drop it
    if max_len_delete:
        # drop it
        of.extend([1e99])
        on.extend([on[-1]])
    else:
        # include it
        of.extend([ind2[-1]])
    #
    pick = []
    while on[-1] > of[0]:
        while on[0] <= of[0]:
            on.popleft()
        while of[0] < on[0]:
            of.popleft()
        if of[0] - on[0] > max_len:
            if max_len_delete:
                on.popleft()
                continue
            of.appendleft(on[0] + max_len)
        pick.append([on[0], of[0]])
    return np.array(pick, dtype=np.int64)

