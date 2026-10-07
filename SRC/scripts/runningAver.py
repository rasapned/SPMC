import sys
import numpy as np
import math

case_name = sys.argv[1]
case_base_name = case_name.split('.')[0]
iter_mult = int(sys.argv[2])

# initialise global arrays
z_it, FEcount_it, O2count_it, O1count_it, H2Ocount_it, H2count_it, H1count_it, OHcount_it, EVcount_it, total_it, ratio_final_it = np.zeros((11, 0), dtype=float)

time_tot, x_tot, Tp_tot, Dp_tot, Tbf_tot = np.zeros((5, 0), dtype=float)

# iterate over all MC runs
for j in range(iter_mult):
    print('Reading and averaging ' + case_name + ' ' + str(j))
    
    # initialise current MC run arrays
    FEcount, O2count, O1count, H2Ocount, H2count, H1count, OHcount, EVcount, total, ratio_final = np.zeros((10, 0), dtype=float)
    time_it, x_it, Tp_it, Dp_it, Tbf_it, delta_Tp, delta_Dp = np.zeros((7, 0), dtype=float)
    z=[]
    
    # read data from file
    file_name = 'RESULTS/' + case_base_name + '_results/' + case_base_name + '-SPMC-' + str(j + 1) + '.csv'
    time, DFB, Tp, Dp, Tbf = np.genfromtxt(file_name, delimiter=',', comments='#', usecols=(0, 1, 3, 4, 14), unpack=True)
    
    # data sampling interval
    deltaT = 1e-8

    # get the position of simulation start (has to be the same for all runs)
    t_cur = math.ceil(time[0] * 100000) / 100000
    if j != 0:
        if t_cur != time_tot[0]:
            print('fuck you')
            exit()

    # interpolate read data at the positions defined by a stationary sampling interval
    time_it = np.arange(t_cur, time[-1], deltaT)
    x_it = np.interp(time_it, time, DFB)
    Tp_it = np.interp(time_it, time, Tp)
    Dp_it = np.interp(time_it, time, Dp)
    Tbf_it = np.interp(time_it, time, Tbf)
    
    # for the first run, just copy read data to global arrays
    if j == 0:
        time_tot = time_it
        x_tot = x_it
        Tp_tot = Tp_it
        Dp_tot = Dp_it
        Tbf_tot = Tbf_it
        # start counting active runs at each DFB
        divisor = np.ones((len(time_tot)))
        sigma_Tp = np.zeros((len(time_tot)))
        sigma_Dp = np.zeros((len(time_tot)))
    
    # for any subsequent runs, adjust the array length (MC sim. can be longer or shorter) and add values to the global array
    else:
        len_new = len(time_it)
        len_old = len(time_tot)
        if len_new > len_old:
            time_tot = np.concatenate((time_tot,time_it[len_old:]))
            x_tot = np.concatenate((x_tot,x_it[len_old:]))
            divisor += 1
            delta_Tp = Tp_it[:len_old] - Tp_tot
            delta_Dp = Dp_it[:len_old] - Dp_tot
            delta_Tbf = Tbf_it[:len_old] - Tbf_tot
            Tp_tot += delta_Tp / divisor
            Dp_tot += delta_Dp / divisor
            Tbf_tot += delta_Tbf / divisor
            sigma_Tp += delta_Tp * (Tp_it[:len_old] - Tp_tot)
            sigma_Dp += delta_Dp * (Dp_it[:len_old] - Dp_tot)
            Tp_tot = np.concatenate((Tp_tot,Tp_it[len_old:]))
            Dp_tot = np.concatenate((Dp_tot,Dp_it[len_old:]))
            Tbf_tot = np.concatenate((Tbf_tot,Tbf_it[len_old:]))
            divisor = np.concatenate((divisor, np.ones((len_new - len_old))))
            sigma_Tp = np.concatenate((sigma_Tp, np.zeros((len_new - len_old))))
            sigma_Dp = np.concatenate((sigma_Dp, np.zeros((len_new - len_old))))

        else:
            
            divisor[:len_new] += 1
            delta_Tp = Tp_it - Tp_tot[:len_new]
            delta_Dp = Dp_it - Dp_tot[:len_new]
            delta_Tbf = Tbf_it - Tbf_tot[:len_new]
            Tp_tot[:len_new] += delta_Tp / divisor[:len_new]
            Dp_tot[:len_new] += delta_Dp / divisor[:len_new]
            Tbf_tot[:len_new] += delta_Tbf / divisor[:len_new]
            sigma_Tp[:len_new] += delta_Tp * (Tp_it - Tp_tot[:len_new])
            sigma_Dp[:len_new] += delta_Dp * (Dp_it - Dp_tot[:len_new])

    # continue the same procedure for reactive collision counts
    with open(file_name, 'r') as f:
        lines = f.readlines()
        SOIlist = []
        ratio = []
        for l in lines:
            if len(l) < 100:
                print(len(SOIlist))
            if l.split(',')[13] == 'el':
                continue
            s = l.split(',')[12]
            SOIlist.append(s.strip())
            if len(SOIlist) > 1:
                ratio.append(float(l.split(',')[6]))

    SOIlist = np.array(SOIlist)
    
    averaging_step = 1e-5
    x_thresh = DFB[0] 
    target = np.arange(x_thresh, DFB[-1], averaging_step)
    browsed = np.asarray(DFB)
    indx_arr = []
    for t in target:
        indx_arr.append((np.abs(browsed - t)).argmin())
    diff_arr = np.asarray(time)[indx_arr]
    timestep = np.diff(diff_arr)
    
    for i in range(len(indx_arr) - 1):
        z.append(DFB[int((indx_arr[i]+indx_arr[i+1])/2)])
        FEcount = np.append(FEcount, np.count_nonzero(SOIlist[indx_arr[i]: indx_arr[i+1]] == 'Fe')/timestep[i])
        O2count = np.append(O2count, np.count_nonzero(SOIlist[indx_arr[i]: indx_arr[i+1]] == 'O2')/timestep[i])
        O1count = np.append(O1count, np.count_nonzero(SOIlist[indx_arr[i]: indx_arr[i+1]] == 'O')/timestep[i]) 
        H2Ocount = np.append(H2Ocount, np.count_nonzero(SOIlist[indx_arr[i]: indx_arr[i+1]] == 'H2O')/timestep[i])
        H2count = np.append(H2count, np.count_nonzero(SOIlist[indx_arr[i]: indx_arr[i+1]] == 'H2')/timestep[i])
        H1count = np.append(H1count, np.count_nonzero(SOIlist[indx_arr[i]: indx_arr[i+1]] == 'H')/timestep[i])
        OHcount = np.append(OHcount, np.count_nonzero(SOIlist[indx_arr[i]: indx_arr[i+1]] == 'OH')/timestep[i]) 
        EVcount = np.append(EVcount, np.count_nonzero(SOIlist[indx_arr[i]: indx_arr[i+1]] == '---')/timestep[i]) 
        total = np.append(total, FEcount[-1] + O2count[-1] + O1count[-1] + H2Ocount[-1] + H2count[-1] + H1count[-1] + OHcount[-1])
        ratio_sum = np.sum(ratio[indx_arr[i]: indx_arr[i + 1]])
        ratio_final = np.append(ratio_final, ratio_sum / (indx_arr[i + 1] - indx_arr[i]))

    if j == 0:
        z_it = z
        FEcount_it = FEcount
        O2count_it = O2count
        O1count_it = O1count
        H2Ocount_it = H2Ocount
        H2count_it = H2count
        H1count_it = H1count
        OHcount_it = OHcount
        EVcount_it = EVcount
        total_it = total
        ratio_final_it = ratio_final
        divisor_bis = np.ones((len(z_it)))
    else:
        len_new = len(z)
        len_old = len(z_it)

        if len_new > len_old:
            z_it = np.concatenate((z_it,z[len_old:]))
            FEcount_it += FEcount[:len_old]
            O2count_it += O2count[:len_old]
            O1count_it += O1count[:len_old]
            H2Ocount_it += H2Ocount[:len_old]
            H2count_it += H2count[:len_old]
            H1count_it += H1count[:len_old]
            OHcount_it += OHcount[:len_old]
            EVcount_it += EVcount[:len_old]
            total_it += total[:len_old]
            ratio_final_it += ratio_final[:len_old]
            divisor_bis += 1

            FEcount_it = np.concatenate((FEcount_it,FEcount[len_old:]))
            O2count_it = np.concatenate((O2count_it,O2count[len_old:]))
            O1count_it = np.concatenate((O1count_it,O1count[len_old:]))
            H2Ocount_it = np.concatenate((H2Ocount_it,H2Ocount[len_old:]))
            H2count_it = np.concatenate((H2count_it,H2count[len_old:]))
            H1count_it = np.concatenate((H1count_it,H1count[len_old:]))
            OHcount_it = np.concatenate((OHcount_it,OHcount[len_old:]))
            EVcount_it = np.concatenate((EVcount_it,EVcount[len_old:]))
            total_it = np.concatenate((total_it,total[len_old:]))
            ratio_final_it = np.concatenate((ratio_final_it,ratio_final[len_old:]))
            divisor_bis = np.concatenate((divisor_bis, np.ones((len_new - len_old))))

        else:
            FEcount_it[:len_new] += FEcount
            O2count_it[:len_new] += O2count
            O1count_it[:len_new] += O1count
            H2Ocount_it[:len_new] += H2Ocount
            H2count_it[:len_new] += H2count
            H1count_it[:len_new] += H1count
            OHcount_it[:len_new] += OHcount
            EVcount_it[:len_new] += EVcount
            total_it[:len_new] += total
            ratio_final_it[:len_new] += ratio_final
            divisor_bis[:len_new] += 1

FEcount_it /= divisor_bis
O2count_it /= divisor_bis
O1count_it /= divisor_bis
H2Ocount_it /= divisor_bis
H2count_it /= divisor_bis
H1count_it /= divisor_bis
OHcount_it /= divisor_bis
EVcount_it /= divisor_bis
total_it /= divisor_bis
ratio_final_it /= divisor_bis


mask = divisor > 1
sigma_Tp_final = np.zeros_like(sigma_Tp)
sigma_Dp_final = np.zeros_like(sigma_Dp)

sigma_Tp_final[mask] = np.sqrt(sigma_Tp[mask] / (divisor[mask] - 1))
sigma_Dp_final[mask] = np.sqrt(sigma_Dp[mask] / (divisor[mask] - 1))
sigma_Tp_final = sigma_Tp / Tp_tot * 100
sigma_Dp_final = sigma_Dp / Dp_tot * 100

countfile_name = 'RESULTS/' + case_base_name + '_results/' + case_base_name + '-count.csv'
with open(countfile_name, 'w') as countfile:
    for i in range(len(z_it)):
        count_str = f"{z_it[i]:1.8e}, {FEcount_it[i]:1.8e}, {O2count_it[i]:1.8e}, {O1count_it[i]:1.8e}, {H2Ocount_it[i]:1.8e}, {H2count_it[i]:1.8e}, {H1count_it[i]:1.8e}, {OHcount_it[i]:1.8e}, {EVcount_it[i]:1.8e}, {total_it[i]:1.8e}, {ratio_final_it[i]:1.8e}\n"
        countfile.write(count_str)

resultfile_name = 'RESULTS/' + case_base_name + '_results/' + case_base_name + '-average.csv'
with open(resultfile_name, 'w') as resultfile:
    for i in range(len(time_tot)):
        result_str = f"{time_tot[i]:1.8e}, {x_tot[i]:1.8e}, {Tp_tot[i]:1.8e}, {Dp_tot[i]:1.8e}, {sigma_Tp[i]:1.8e}, {sigma_Dp[i]:1.8e}, {Tbf_tot[i]:1.8e}\n"
        resultfile.write(result_str)

print('Done!')

