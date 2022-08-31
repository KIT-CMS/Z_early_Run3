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
    print("Group of Runs: ", run_group_list)
    return run_group_list
