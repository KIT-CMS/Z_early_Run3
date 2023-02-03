#function that groups runs by luminosity set by the lumi_limit_for_grouping
def run_lumi_selection(prev_run_number, run_list, run_lumi, lumi_limit_for_grouping):
    sum_lumi = 0
    run_group_list = list()
    i = 0
    if len(run_list) == 1:
        run_group_list.append(run_list[0])
        return run_group_list
    while i < len(run_list):
        if float(prev_run_number) != float(run_list[0]):
            if i == 0:
                while float(run_list[i]) < float(prev_run_number):
                    i+=1
                i+=1
        if len(run_group_list) == 0 and float (run_lumi[run_list[i]]) > lumi_limit_for_grouping:
                run_group_list = list()
                run_group_list.append(run_list[i])
                return run_group_list
        temp_var = run_lumi[run_list[i]]
        sum_lumi+=float(temp_var)
        if sum_lumi >= lumi_limit_for_grouping:
            if float (run_lumi[run_list[i]]) > lumi_limit_for_grouping:
                break
            else:
                run_group_list.append(run_list[i])
                break
        else:
            run_group_list.append(run_list[i])
        i+=1
    # print("Group of Runs: ", run_group_list, sum_lumi)
    return run_group_list

def get_subera_boundaries():
    subera_boundaries = {
        "Run2022": [
            (355862, 357900),
        ],
        "C": [
            (355862, 357482),
        ],
        "D1": [
            (357538, 357732),
            # (357483, 357733),
        ],
        "D2": [
            (357734, 357900),
        ],
    }

    # 355862-356446 660.0894283320001
    # 356523-357482 4184.218497299999

    # 355862-357482 4844.307925632
    # 357538-357732 919.2951682019999
    # 357734-357900 1781.2130120119998

    # (355862, 356446),  # early C
    # (356446+1, 357482), # later C

    # (355862, 357482),  # C
    # (357482+1, 357733), # D v1
    # (357733+1, 357900), # D v2

    return subera_boundaries

def run_by_era(run_list, subera):
    subera_boundaries = get_subera_boundaries()
    out = []
    for run_i, run_f in subera_boundaries[subera]:
        run_group_list = []
        for run in run_list:
            if int(run) >= int(run_i) and int(run) <= int(run_f):
                run_group_list.append(run)
        out.append(run_group_list)

    return out
