from pickletools import float8
from ntuple_processor.utils import Selection

#if grouping runs by luminosity this returns the group
def run_group_selection(initial_run, final_run):
    return Selection(name=initial_run + "-" + final_run, cuts=[("run >= {} && run <= {}".format(initial_run, final_run), "run_{}_{}".format(initial_run, final_run))])

def run_selection(_run):
    return Selection(name=_run, cuts=[("run == {}".format(_run), _run)])