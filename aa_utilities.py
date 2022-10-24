#%%
# This file contains utility functions for the calculations in
# Matheson/Malhotra/Keane 'A von Mises-Fisher distribution for the Plutinos'.
#%%
def add_category_lists_sv20classifer(int_dict,prediction,Detached_list,Resonant_list,Scattering_list,Classical_list):
# Takes in a dictionary of class names.
# Takes in a prediction (probabilities that an object belongs to each class).
# Takes in lists of probabilities assigned to membership of each class, for each prior object.
# Sorts out the probabilities and assigns them to the correct list,
# ie, Classical probability to the Classical list, Resonant probability to the Resonant list, etc.
    if int_dict[0] == 'Detached':
        Detached_list.append(prediction[0][0])
    elif int_dict[0] == 'Resonant':
        Resonant_list.append(prediction[0][0])
    elif int_dict[0] == 'Scattering':
        Scattering_list.append(prediction[0][0])
    elif int_dict[0] == 'Classical':
        Classical_list.append(prediction[0][0])
    else:
        raise Exception('Unknown format in int_dict.',int_dict)
    if int_dict[1] == 'Detached':
        Detached_list.append(prediction[0][1])
    elif int_dict[1] == 'Resonant':
        Resonant_list.append(prediction[0][1])
    elif int_dict[1] == 'Scattering':
        Scattering_list.append(prediction[0][1])
    elif int_dict[1] == 'Classical':
        Classical_list.append(prediction[0][1])
    else:
        raise Exception('Unknown format in int_dict.',int_dict)
    if int_dict[2] == 'Detached':
        Detached_list.append(prediction[0][2])
    elif int_dict[2] == 'Resonant':
        Resonant_list.append(prediction[0][2])
    elif int_dict[2] == 'Scattering':
        Scattering_list.append(prediction[0][2])
    elif int_dict[2] == 'Classical':
        Classical_list.append(prediction[0][2])
    else:
        raise Exception('Unknown format in int_dict.',int_dict)
    if int_dict[3] == 'Detached':
        Detached_list.append(prediction[0][3])
    elif int_dict[3] == 'Resonant':
        Resonant_list.append(prediction[0][3])
    elif int_dict[3] == 'Scattering':
        Scattering_list.append(prediction[0][3])
    elif int_dict[3] == 'Classical':
        Classical_list.append(prediction[0][3])
    else:
        raise Exception('Unknown format in int_dict.',int_dict)
    return Detached_list,Resonant_list,Scattering_list,Classical_list
#%%
def circular_hist(ax, x, bins=16, density=True, offset=0, gaps=True):
    """
    https://stackoverflow.com/questions/22562364/circular-polar-histogram-in-python
    Produce a circular histogram of angles on ax.
    Parameters
    ----------
    ax : matplotlib.axes._subplots.PolarAxesSubplot
        axis instance created with subplot_kw=dict(projection='polar').
    x : array
        Angles to plot, expected in units of radians.
    bins : int, optional
        Defines the number of equal-width bins in the range. The default is 16.
    density : bool, optional
        If True plot frequency proportional to area. If False plot frequency
        proportional to radius. The default is True.
    offset : float, optional
        Sets the offset for the location of the 0 direction in units of
        radians. The default is 0.
    gaps : bool, optional
        Whether to allow gaps between bins. When gaps = False the bins are
        forced to partition the entire [-pi, pi] range. The default is True.
    Returns
    -------
    n : array or list of arrays
        The number of values in each bin.
    bins : array
        The edges of the bins.
    patches : `.BarContainer` or list of a single `.Polygon`
        Container of individual artists used to create the histogram
        or list of such containers if there are multiple input datasets.
    """
    import numpy as np
    # Wrap angles to [-pi, pi)
    x = (x+np.pi) % (2*np.pi) - np.pi
    # Force bins to partition entire circle
    if not gaps:
        bins = np.linspace(-np.pi, np.pi, num=bins+1)
    # Bin data and record counts
    n, bins = np.histogram(x, bins=bins)
    # Compute width of each bin
    widths = np.diff(bins)
    # By default plot frequency proportional to area
    if density:
        # Area to assign each bin
        area = n / x.size
        # Calculate corresponding bin radius
        radius = (area/np.pi) ** .5
    # Otherwise plot frequency proportional to radius
    else:
        radius = n
    # Plot data on ax
    patches = ax.bar(bins[:-1], radius, zorder=1, align='edge', width=widths,
                     edgecolor='C0', fill=False, linewidth=1)
    # Set the direction of the zero angle
    ax.set_theta_offset(offset)
    # Remove ylabels for area plots (they are mostly obstructive)
    if density:
        ax.set_yticks([])
    return n, bins, patches
#%% check if object is resonant over 0-tmax time period
# Resonance identification algorithm described in
# Yu/Murray-Clay/Volk 2018, AJ 156:33, https://doi.org/10.3847/1538-3881/aac6cd
def check_resonance(orbels_dict,settings):
    # settings = {'tmax_yrs':tmax_yrs,'tstep_yrs':tstep_yrs,'min_stick_length_yrs':\
    # min_stick_length_yrs,'min_stick_length':min_stick_length,'bigmax':bigmax,\
    # 'smallmax':smallmax,'ordermax':ordermax,'smatol':smatol,'shift_amount':\
    # shift_amount,'mmr_tolerance':mmr_tolerance,'libration_threshold':libration_threshold,\
    # 'points_to_average':points_to_average,'averaging_proportion':averaging_proportion,\
    # 'long_or_small_threshold_pt_count':long_or_small_threshold_pt_count,\
    # 'small_division_count':small_division_count,'long_division_count':\
    # long_division_count,'length_threshold':length_threshold,'scattering_tolerance_delta_P':\
    # scattering_tolerance_delta_P,'scattering_tolerance_delta_a':scattering_tolerance_delta_a,\
    # 'tjmax':tjmax,'qmax':qmax,'a_oort':a_oort,'e_detached':e_detached,'a_outer':a_outer,\
    # 'a_inner':a_inner,'infile':infile,'outfile':outfile,'settingsfile':settingsfile,\
    # 'starting_sim_file':starting_sim_file,'integrator':integrator,'JD':JD}

    # orbels_dict = {'t':t_yrs_array,'a':a_array,'a_Jupiter':a_Jupiter_array,\
    # 'a_Neptune':a_Neptune_array,'tj':tj_array,'q':q_array,'e':e_array,\
    # 'inc':inc_array,'w':w_array,'W':W_array,'M':M_array,'l':l_array,\
    # 'l_Neptune':l_Neptune_array,'pomega':pomega_array,'pomega_Neptune':\
    # pomega_Neptune_array,'f':f_array,'f_Neptune':f_Neptune_array,\
    # 'P':P_yrs_array,'P_Neptune':P_yrs_Neptune_array,'theta':theta_array,\
    # 'theta_Neptune':theta_Neptune_array}
    import numpy as np
    bigmax = settings['bigmax'][0]
    smallmax = settings['smallmax'][0]
    ordermax = settings['ordermax'][0]
    min_stick_length = settings['min_stick_length'][0]
    shift_amount = settings['shift_amount'][0]
    scattering_tolerance_delta_a = settings['scattering_tolerance_delta_a'][0]
    scattering_tolerance_delta_P = settings['scattering_tolerance_delta_P'][0]
    mmr_tolerance = settings['mmr_tolerance'][0]
    points_to_average = settings['points_to_average'][0]
    averaging_proportion = settings['averaging_proportion'][0]
    libration_threshold = settings['libration_threshold'][0]
    long_or_small_threshold_pt_count = settings['long_or_small_threshold_pt_count'][0]
    long_division_count = settings['long_division_count'][0]
    small_division_count = settings['small_division_count'][0]
    length_threshold = settings['length_threshold'][0]
    a_list_Neptune = orbels_dict['a_Neptune']
    l_pluto_array = orbels_dict['l']
    l_Neptune_array = orbels_dict['l_Neptune']
    pomega_pluto_array = orbels_dict['pomega']
    a_list_1 = orbels_dict['a']
    amax = 150
    a_Neptune = a_list_Neptune[0]
    amin = a_Neptune
    df_mmrs = make_mmrs(amin,amax,a_Neptune,bigmax,smallmax,ordermax)
    Npts = len(l_pluto_array)
    Nmmr = df_mmrs.shape[0]
    startpts_list = np.arange(0,Npts-min_stick_length,shift_amount)
    stoppts_list = startpts_list + min_stick_length
    if stoppts_list[-1] < Npts:
        stoppts_list[-1] = Npts # add a couple extra points on the end
    Nwindow = len(startpts_list)
    resonant_window_list = []
    resonant_big_list = []
    resonant_small_list = []
    resonant_diff_list = []
    for iwin in range(Nwindow):
        any_resonance = False
        start_pt = startpts_list[iwin]
        stop_pt = stoppts_list[iwin]
        # do stuff in this sliding window
        a_list_2 = a_list_1[start_pt:stop_pt]
        a_list_Neptune_2 = a_list_Neptune[start_pt:stop_pt]
        l_list_2 = l_pluto_array[start_pt:stop_pt]
        l_list_Neptune_2 = l_Neptune_array[start_pt:stop_pt]
        pomega_list_2 = pomega_pluto_array[start_pt:stop_pt]
        # first check that particle isn't actively scattering
        a_list_2 = np.array(a_list_2)
        a_list_Neptune_2 = np.array(a_list_Neptune_2)
        l_list_2 = np.array(l_list_2)
        l_list_Neptune_2 = np.array(l_list_Neptune_2)
        pomega_list_2 = np.array(pomega_list_2)
        if np.min(a_list_2) > 0:
            period_ratio_list_2 = (a_list_2/a_list_Neptune_2)**(3/2)
            pmax = np.max(period_ratio_list_2)
            pmin = np.min(period_ratio_list_2)
            delta_p = (pmax-pmin)/pmin
            delta_a = np.max(a_list_2) - np.min(a_list_2)
        else:
            delta_p = 99999 # arbitrarily large number
        if (delta_p > scattering_tolerance_delta_P) or \
           (delta_a > scattering_tolerance_delta_a):
            2-2 # move on to the next window
        else: # do stuff
            2-2
            # check which resonances to look for
            pmean = np.mean(period_ratio_list_2)
            for immr in range(Nmmr):
                pmmr = df_mmrs['period_ratio'][immr]
                # print(iwin+1,Nwindow,immr+1,Nmmr,pmmr)
                if abs(pmmr-pmean) <= mmr_tolerance:
                    # calculate resonant angle for this mmr in this window
                    big = df_mmrs['big'][immr]
                    small = df_mmrs['small'][immr]
                    theta_list_1 = big * l_list_2 - small * l_list_Neptune_2 - (big-small) * pomega_list_2
                    theta_list_1 = np.mod(theta_list_1,360)
                    # check for libration
                    theta_sorted = np.sort(theta_list_1)
                    smallest_ten = theta_sorted[0:points_to_average]
                    largest_ten = theta_sorted[-points_to_average:]
                    smallest_mean = np.mean(smallest_ten)
                    largest_mean = np.mean(largest_ten)
                    libration_amplitude = (largest_mean-smallest_mean)/2
                    # print(libration_amplitude)
                    if libration_amplitude > libration_threshold:
                        librates_around_180 = False
                    else:
                        librates_around_180 = True
                        libration_amplitude_180 = libration_amplitude
                    theta_list_2 = theta_list_1 - 360 * (theta_list_1>180) # shift from 0,360 to -180,180
                    theta_sorted = np.sort(theta_list_2)
                    smallest_ten = theta_sorted[0:points_to_average]
                    largest_ten = theta_sorted[-points_to_average:]
                    smallest_mean = np.mean(smallest_ten)
                    largest_mean = np.mean(largest_ten)
                    libration_amplitude = (largest_mean-smallest_mean)/2
                    # print(libration_amplitude)
                    if libration_amplitude > libration_threshold:
                        librates_around_0 = False
                    else:
                        librates_around_0 = True
                        libration_amplitude_0 = libration_amplitude
                    if librates_around_180 or librates_around_0:
                    # if librates_around_180:
                        if librates_around_0:
                            libration_amplitude = libration_amplitude_0
                        if librates_around_180:
                            libration_amplitude = libration_amplitude_180
                        # record which time window and which resonance showed libration
                        if len(resonant_window_list) == 0:
                            resonant_window_list.append(iwin)
                        elif resonant_window_list[-1] != iwin:
                            resonant_window_list.append(iwin)
                        if len(resonant_big_list) == 0:
                            resonant_big_list.append([big])
                        elif len(resonant_big_list) != len(resonant_window_list):
                            resonant_big_list.append([big])
                        else:
                            resonant_big_list[-1].append(big)
                        if len(resonant_small_list) == 0:
                            resonant_small_list.append([small])
                        elif len(resonant_small_list) != len(resonant_window_list):
                            resonant_small_list.append([small])
                        else:
                            resonant_small_list[-1].append(small)
                        diff = big - small
                        if len(resonant_diff_list) == 0:
                            resonant_diff_list.append([diff])
                        elif len(resonant_diff_list) != len(resonant_window_list):
                            resonant_diff_list.append([diff])
                        else:
                            resonant_diff_list[-1].append(diff)
                        any_resonance = True
            if (any_resonance == False):
                2-2
    '''
    identify the unique resonances in the resonance list, then find the
    starts/stops/lengths/centers/etc of each resonance over time
    '''
    unique_big_list = []
    unique_small_list = []
    unique_diff_list = []
    unique_combo_list = []
    Nreswin = len(resonant_window_list)
    Nbigwin = len(resonant_big_list)
    Nsmallwin = len(resonant_small_list)
    Ndiffwin = len(resonant_diff_list)
    listwin = [Nreswin,Nbigwin,Nsmallwin,Ndiffwin]
    if np.std(listwin) != 0.0:
        2-2
    for iwin2 in range(Nreswin):
        big_list_here = resonant_big_list[iwin2]
        small_list_here = resonant_small_list[iwin2]
        diff_list_here = resonant_diff_list[iwin2]
        if (len(big_list_here) != len(small_list_here)) or (len(big_list_here) != len(diff_list_here)):
            2-2
        Nres_here = len(big_list_here)
        for ires in range(Nres_here):
            big_here = big_list_here[ires]
            small_here = small_list_here[ires]
            diff_here = diff_list_here[ires]
            combo = str(big_here) + ' to ' + str(small_here)
            if combo not in unique_combo_list:
                unique_combo_list.append(combo)
                unique_big_list.append(big_here)
                unique_small_list.append(small_here)
                unique_diff_list.append(diff_here)
    Nunique = len(unique_combo_list)
    unique_sticks_windows_lists = []
    unique_continuous_sticks_points_start_lists = []
    unique_continuous_sticks_points_end_lists = []
    unique_continuous_sticks_proportionatelength_lists = []
    unique_sticks_totalproportionatelength_list = []
    unique_continuous_sticks_amplitude_lists = []
    unique_continuous_sticks_center_lists = []
    for ires2 in range(Nunique):
        big = unique_big_list[ires2]
        small = unique_small_list[ires2]
        sticks_windows_list = []
        continuous_sticks_windows_start_list = []
        continuous_sticks_windows_end_list = []
        continuous_sticks_amplitude_list = []
        continuous_sticks_proportionatelength_list = []
        continuous_sticks_center_list = []
        for iwin3 in range(Nreswin):
            big_list_here = resonant_big_list[iwin3]
            small_list_here = resonant_small_list[iwin3]
            if (big in big_list_here) and (small in small_list_here):
                sticks_windows_list.append(iwin3)
        unique_sticks_windows_lists.append(sticks_windows_list)
        Nstickwin = len(sticks_windows_list)
        # print('\n')
        for iwin4 in range(Nstickwin):
            stickwin = sticks_windows_list[iwin4]
            if len(continuous_sticks_windows_start_list) == 0:
                continuous_sticks_windows_start_list.append(stickwin)
            if (iwin4 != 0) and (stickwin>(sticks_windows_list[iwin4-1]+1)):
                continuous_sticks_windows_start_list.append(stickwin)
            if (iwin4 == Nstickwin-1):
                continuous_sticks_windows_end_list.append(stickwin)
            elif (stickwin<(sticks_windows_list[iwin4+1]-1)):
                continuous_sticks_windows_end_list.append(stickwin)
        startpts_list_2 = []
        stoppts_list_2 = []
        Nstickwin2 =  len(continuous_sticks_windows_start_list)
        for iwin6 in range(Nstickwin2):
            startpt = startpts_list[continuous_sticks_windows_start_list[iwin6]]
            stoppt = stoppts_list[continuous_sticks_windows_end_list[iwin6]]
            startpts_list_2.append(startpt)
            stoppts_list_2.append(stoppt)
        stoppts_list_3 = []
        stopping_yet = False
        stop_index = -1
        while_ct_1 = 0
        while stopping_yet == False:
            while_ct_1 = while_ct_1 + 1
            stop_index = stop_index + 1
            stoppt = stoppts_list_2[stop_index]
            if stop_index == Nstickwin2-1:
                stopping_yet = True
            else:
                next_startpt = startpts_list_2[stop_index+1]
                if next_startpt <= stoppt:
                    stopping_yet = False
                else:
                    stopping_yet = True
        stoppts_list_3.append(stoppts_list_2[stop_index])
        while_ct_2 = 0
        while stop_index != Nstickwin2-1:
            while_ct_2 = while_ct_2 + 1
            stopping_yet = False
            inner_while_ct_2 = 0
            while stopping_yet == False:
                inner_while_ct_2 = inner_while_ct_2 + 1
                stop_index = stop_index + 1
                stoppt = stoppts_list_2[stop_index]
                if stop_index == Nstickwin2-1:
                    stopping_yet = True
                else:
                    next_startpt = startpts_list_2[stop_index+1]
                    if next_startpt <= stoppt:
                        stopping_yet = False
                    else:
                        stopping_yet = True
            stoppts_list_3.append(stoppts_list_2[stop_index])
        startpts_list_3 = [startpts_list_2[0]]
        start_index = 0
        while_ct_3 = 0
        while len(startpts_list_3) < len(stoppts_list_3):
            while_ct_3 = while_ct_3 + 1
            start_index = start_index + 1
            if start_index < len(startpts_list_2):
                if startpts_list_2[start_index] > stoppts_list_2[start_index-1]:
                    startpts_list_3.append(startpts_list_2[start_index])
            else:
                break
        unique_continuous_sticks_points_start_lists.append(startpts_list_3)
        unique_continuous_sticks_points_end_lists.append(stoppts_list_3)
        Nconwin = len(startpts_list_3)
        for iwin5 in range(Nconwin):
            startpt = startpts_list_3[iwin5]
            stoppt = stoppts_list_3[iwin5]
            Npts_here = stoppt-startpt
            proportionatelength = Npts_here/Npts
            continuous_sticks_proportionatelength_list.append(proportionatelength)
        unique_continuous_sticks_proportionatelength_lists.append(\
                continuous_sticks_proportionatelength_list)
        unique_sticks_totalproportionatelength_list.append(\
                np.sum(continuous_sticks_proportionatelength_list))
        # # # recompute amplitudes and centers of continuous sticks
        # # # this should eliminate spurious (circulating) angles
        Nstart = len(startpts_list_3)
        Nstop = len(stoppts_list_3)
        median_amplitude_list = []
        median_center_list = []
        if Nstart != Nstop:
            2-2
        for istart in range(Nstart):
            startpt = startpts_list_3[istart]
            stoppt = stoppts_list_3[istart]
            l_pluto_here = l_pluto_array[startpt:stoppt]
            l_Neptune_here = l_Neptune_array[startpt:stoppt]
            pomega_pluto_here = pomega_pluto_array[startpt:stoppt]
            theta_list_here = big * l_pluto_here - small * l_Neptune_here - (big-small) * pomega_pluto_here
            theta_list_here = np.mod(theta_list_here,360)
            Ntheta = len(theta_list_here)
            if Ntheta >= long_or_small_threshold_pt_count:
                division_count = long_division_count
            else:
                division_count = small_division_count
            libration_amplitude_list = []
            libration_center_list = []
            Npts_per = int(np.floor(Ntheta/division_count))
            for idiv in range(division_count):
                start_pt = Npts_per*idiv
                stop_pt = start_pt + Npts_per
                if abs(stop_pt-Ntheta) < Npts_per:
                    stop_pt = Ntheta
                theta_here = theta_list_here[start_pt:stop_pt]
                theta_sorted = np.sort(theta_here)
                Npts_here = len(theta_sorted)
                points_to_average = int(np.ceil(Npts_here*averaging_proportion))
                smallest_few = theta_sorted[0:points_to_average]
                largest_few = theta_sorted[-points_to_average:]
                smallest_mean = np.mean(smallest_few)
                largest_mean = np.mean(largest_few)
                libration_amplitude = (largest_mean-smallest_mean)/2
                libration_amplitude_list.append(libration_amplitude)
                libration_center = np.mean([largest_mean,smallest_mean])
                libration_center_list.append(libration_center)
            median_amplitude = np.median(libration_amplitude_list)
            median_center = np.median(libration_center_list)
            median_amplitude_list.append(median_amplitude)
            median_center_list.append(median_center)
        unique_continuous_sticks_amplitude_lists.append(median_amplitude_list)
        unique_continuous_sticks_center_lists.append(median_center_list)
    unique_sticks_windows_lists_alt = []
    for ires4 in range(Nunique):
        sticks_windows_list = unique_sticks_windows_lists[ires4]
        sticks_windows_lists = []
        Nswl = len(sticks_windows_list)
        iswl = 0
        while iswl < Nswl:
            win = sticks_windows_list[iswl]
            if len(sticks_windows_lists) == 0:
                sticks_windows_lists.append([win])
            elif win == sticks_windows_lists[-1][-1]+1:
                sticks_windows_lists[-1].append(win)
            else:
                sticks_windows_lists.append([win])
            iswl = iswl + 1
        unique_sticks_windows_lists_alt.append(sticks_windows_lists)
    unique_combo_list_2 = []
    unique_big_list_2 = []
    unique_small_list_2 = []
    unique_diff_list_2 = []
    unique_sticks_windows_lists_2 = []
    unique_continuous_sticks_points_start_lists_2 = []
    unique_continuous_sticks_points_end_lists_2 = []
    unique_continuous_sticks_proportionatelength_lists_2 = []
    unique_sticks_totalproportionatelength_list_2 = []
    unique_continuous_sticks_amplitude_lists_2 = []
    unique_continuous_sticks_center_lists_2 = []
    for ires3 in range(Nunique):
        combo = unique_combo_list[ires3]
        big = unique_big_list[ires3]
        small = unique_small_list[ires3]
        diff = unique_diff_list[ires3]
        windows_lists = unique_sticks_windows_lists_alt[ires3]
        continuous_sticks_points_start_list = unique_continuous_sticks_points_start_lists[ires3]
        continuous_sticks_points_end_list = unique_continuous_sticks_points_end_lists[ires3]
        continuous_sticks_proportionatelength_list = unique_continuous_sticks_proportionatelength_lists[ires3]
        continuous_sticks_amplitude_list = unique_continuous_sticks_amplitude_lists[ires3]
        continuous_sticks_center_list = unique_continuous_sticks_center_lists[ires3]
        Nsticks = len(continuous_sticks_points_start_list)
        windows_lists_2 = []
        start_list_2 = []
        stop_list_2 = []
        length_list_2 = []
        amplitude_list_2 = []
        center_list_2 = []
        for istick in range(Nsticks):
            windows = windows_lists[istick]
            startpt = continuous_sticks_points_start_list[istick]
            stoppt = continuous_sticks_points_end_list[istick]
            length = continuous_sticks_proportionatelength_list[istick]
            amplitude = continuous_sticks_amplitude_list[istick]
            center = continuous_sticks_center_list[istick]
            if amplitude <= libration_threshold:
                windows_lists_2.append(windows)
                start_list_2.append(startpt)
                stop_list_2.append(stoppt)
                length_list_2.append(length)
                amplitude_list_2.append(amplitude)
                center_list_2.append(center)
        Nconwin2 = len(start_list_2)
        if Nconwin2 > 0:
            unique_combo_list_2.append(combo)
            unique_big_list_2.append(big)
            unique_small_list_2.append(small)
            unique_diff_list_2.append(diff)
            unique_sticks_windows_lists_2.append(windows_lists_2)
            unique_continuous_sticks_points_start_lists_2.append(start_list_2)
            unique_continuous_sticks_points_end_lists_2.append(stop_list_2)
            unique_continuous_sticks_proportionatelength_lists_2.append(length_list_2)
            unique_sticks_totalproportionatelength_list_2.append(np.sum(length_list_2))
            unique_continuous_sticks_amplitude_lists_2.append(amplitude_list_2)
            unique_continuous_sticks_center_lists_2.append(center_list_2)
            theta_list_plot = big * l_pluto_array - small * l_Neptune_array - (big-small) * pomega_pluto_array
            theta_list_plot = np.mod(theta_list_plot,360)
    # # identify resonance with longest sticking time, if any
    Nunique2 = len(unique_combo_list_2)
    if Nunique2 == 0: # not resonant
        resonant = False
        big = 'NA'
        small = 'NA'
    else: # in one or more mean motion resonances for some amount of time
        lenlist = unique_sticks_totalproportionatelength_list_2
        lenindex = lenlist.index(np.max(lenlist))
        lenmax = lenlist[lenindex]
        big = unique_big_list_2[lenindex]
        small = unique_small_list_2[lenindex]
        big = str(big)
        small = str(small)
        if lenmax >= length_threshold:
            resonant = True
        else:
            resonant = False
            big = 'NA'
            small = 'NA'
    return resonant,big,small
#%%
def datestr(datevec,flag):
# Takes in a vector of the format [yyyy,mm,dd,hh,mm,ss].
# No leading zeros allowed in the elements of th input vector.
# Returns a string of the format 'yyyymmdd' for flag 'short'.
# Returns a string of the format 'yyyy_mm_dd_hh_mm_ss' for flag 'long'.
    if flag == 'short':
        yyyystr = str(datevec[0])
        mmstr = str(datevec[1])
        ddstr = str(datevec[2])
        if len(mmstr) == 1:
            mmstr = '0' + mmstr
        if len(ddstr) == 1:
            ddstr = '0' + ddstr
        datestr = yyyystr + mmstr + ddstr
    else:
        if flag == 'long':
            datestr = ''
            for i in range(len(datevec)-1):
                thisstr = str(datevec[i])
                if len(thisstr) == 1:
                    thisstr = '0' + thisstr
                if i == 0:
                    datestr = datestr + thisstr
                else:
                    datestr = datestr + '_' + thisstr
            thisstr = str(datevec[-1])
            if datevec[-1] < 10:
                thisstr = '0' + thisstr
            datestr = datestr + '_' + thisstr
    return datestr
#%%
def dat2csv(date):
# All the info we really need from MPCORB is :
# packed designation, semimajor axis, number of oppositions.
    import pandas as pd
    import os.path
    infile = '00_MPCORB_' + date + '.DAT'
    outfile = '00_MPCORB_reduced_' + date + '.csv'
    if os.path.exists(infile) and not os.path.exists(outfile):
        print('In dat2csv, infile exists, making outfile.')
        designation_list = []
        a_list = []
        Nopp_list = []
        file = open(infile,'r')
        lines = file.readlines()
        Nlines = len(lines)
        for iline in range(Nlines):
            if iline > 42:
                line = lines[iline]
                if len(line) > 3: # protect against blank lines
                    designation = line[0:7]
                    designation = designation.lstrip()
                    designation = designation.rstrip()
                    designation_list.append(designation)
                    a = line[92:103]
                    a = a.lstrip()
                    a = a.rstrip()
                    a = float(a)
                    a_list.append(a)
                    Nopp = line[123:126]
                    Nopp = Nopp.lstrip()
                    Nopp = Nopp.rstrip()
                    Nopp = int(Nopp)
                    Nopp_list.append(Nopp)
        this_dict = {'packed_designation':designation_list,'a_au':a_list,'Nopp':Nopp_list}
        df = pd.DataFrame.from_dict(this_dict)
        # Trim the dataframe so we save as small a file as possible.
        # Get gmv08 criteria for which objects we can try to classify.
        gmv08settings = get_gmv08settings()
        min_oppositions = gmv08settings['Nopp_min']
        # Drop objects in the MPC list outside of rough semimajor axis bounds and
        # with too few observed oppositions to be worth classifying.
        df = df[df['a_au']>28] # rough cut
        df = df[df['a_au']<200] # rough cut
        df = df[df['Nopp']>=min_oppositions]
        df.to_csv(outfile,index=False)
    elif os.path.exists(infile) and os.path.exists(outfile):
        print('In dat2csv, infile exists, outfile already exists, doing nothing.')
    elif not os.path.exists(infile) and not os.path.exists(outfile):
        raise Exception('In dat2csv, infile does not exist, outfile does not exist, cannot do anything.')
    elif not os.path.exists(infile) and os.path.exists(outfile):
        print('In dat2csv, infile does not exist, outfile already exists, doing nothing.')
    return
#%%
def generate_clones(JD,datestr,Nclones):
# Read in covariance for each object.
# Generate clones for each object from the covariance.
# For each clone, save heliocentric elements of the object and the planets.
# These will be used later in building a simulation in Rebound for classifying
# the clones. It is easiest to build a simulation from heliocentric elements.
    import numpy as np
    import pandas as pd
    import json
    import urllib
    from astroquery.jplhorizons import Horizons
    # Start with the list of objects left as a rough cut by sbdb_reduce.
    infile = 'sbdb_reduce_output.csv'
    # Start a data frame with the objects left as a rough cut by sbdb_reduce.
    df = pd.read_csv(infile)
    Nobj = df.shape[0]
    center = '500@10' # easiest to build a simulation in heliocentric elements
    for iobj in range(Nobj):
        # initialize lists of heliocentric elements of the object and the planets
        epochP = []
        ePh = [] # eccentricity, Plutino, heliocentric
        qPh_au = [] # perihelion distance, Plutino, heliocentric, au
        tpPh_jd = [] # time of perihelion passage, Plutino, heliocentric, Julian date TDB
        WPh_deg = [] # longitude of ascending node, Plutino, heliocentric, degrees
        wPh_deg = [] # argument of perihelion, Plutino, heliocentric, degrees
        iPh_deg = [] # inclination, Plutino, heliocentric, degrees
        eJh = [] # Jupiter
        qJh_au = []
        tpJh_jd = []
        WJh_deg = []
        wJh_deg = []
        iJh_deg = []
        eSh = [] # Saturn
        qSh_au = []
        tpSh_jd = []
        WSh_deg = []
        wSh_deg = []
        iSh_deg = []
        eUh = [] # Uranus
        qUh_au = []
        tpUh_jd = []
        WUh_deg = []
        wUh_deg = []
        iUh_deg = []
        eNh = [] # Neptune
        qNh_au = []
        tpNh_jd = []
        WNh_deg = []
        wNh_deg = []
        iNh_deg = []
        # Get packed MPC designation of object.
        des = df['packed_designation'][iobj]
        # Contact JPL Solar System Dynamics Small Body Database API to get covariance.
        # Returns heliocentric orbels and covariance
        url = 'https://ssd-api.jpl.nasa.gov/sbdb.api?sstr='+des+'&full-prec=True&cov=mat'
        response = urllib.request.urlopen(url)
        data = json.loads(response.read())
        # Save covariance and other information to file.
        jsonfile = 'ssd_json_' + des + '.json'
        with open(jsonfile,'w') as of:
            json.dump(data,of)
        # Save heliocentric orbital elements of the Plutino and the planets.
        clones_file = 'clones_' + des + '.csv'
        planets_file = 'planets_for_' + des + '.csv'
        if des == 'D4340': # Pluto
            # Because Pluto's covariance is essentially zero, we save heliocentric orbels at
            # the arbitrary Julian date used as an input.
            # Save heliocentric orbital elements of the planets at the arbitrary epoch.
            epochP.append(JD)
            obj = Horizons(id='5',location=center,epochs=JD) # Jupiter barycenter
            el = obj.elements()
            eJh.append(float(el['e']))
            qJh_au.append(float(el['q']))
            tpJh_jd.append(float(el['Tp_jd']))
            WJh_deg.append(np.mod(float(el['Omega']),360))
            wJh_deg.append(np.mod(float(el['w']),360))
            iJh_deg.append(float(el['incl']))
            obj = Horizons(id='6',location=center,epochs=JD) # Saturn barycenter
            el = obj.elements()
            eSh.append(float(el['e']))
            qSh_au.append(float(el['q']))
            tpSh_jd.append(float(el['Tp_jd']))
            WSh_deg.append(np.mod(float(el['Omega']),360))
            wSh_deg.append(np.mod(float(el['w']),360))
            iSh_deg.append(float(el['incl']))
            obj = Horizons(id='7',location=center,epochs=JD) # Uranus barycenter
            el = obj.elements()
            eUh.append(float(el['e']))
            qUh_au.append(float(el['q']))
            tpUh_jd.append(float(el['Tp_jd']))
            WUh_deg.append(np.mod(float(el['Omega']),360))
            wUh_deg.append(np.mod(float(el['w']),360))
            iUh_deg.append(float(el['incl']))
            obj = Horizons(id='8',location=center,epochs=JD) # Neptune barycenter
            el = obj.elements()
            eNh.append(float(el['e']))
            qNh_au.append(float(el['q']))
            tpNh_jd.append(float(el['Tp_jd']))
            WNh_deg.append(np.mod(float(el['Omega']),360))
            wNh_deg.append(np.mod(float(el['w']),360))
            iNh_deg.append(float(el['incl']))
            dictionary = {'epochP_jd':epochP,
                  'eJh':eJh,'qJh_au':qJh_au,'tpJh_jd':tpJh_jd,'WJh_deg':WJh_deg,'wJh_deg':wJh_deg,'iJh_deg':iJh_deg,\
                  'eSh':eSh,'qSh_au':qSh_au,'tpSh_jd':tpSh_jd,'WSh_deg':WSh_deg,'wSh_deg':wSh_deg,'iSh_deg':iSh_deg,\
                  'eUh':eUh,'qUh_au':qUh_au,'tpUh_jd':tpUh_jd,'WUh_deg':WUh_deg,'wUh_deg':wUh_deg,'iUh_deg':iUh_deg,\
                  'eNh':eNh,'qNh_au':qNh_au,'tpNh_jd':tpNh_jd,'WNh_deg':WNh_deg,'wNh_deg':wNh_deg,'iNh_deg':iNh_deg}
            df_out = pd.DataFrame(dictionary)
            df_out.to_csv(planets_file,index=False)
            # The first line of the clones_[des].csv file is always the nominal orbit at epoch.
            obj = Horizons(id='9',location=center,epochs=JD) # Pluto-Charon barycenter
            el = obj.elements()
            ePh.append(float(el['e']))
            qPh_au.append(float(el['q'])) # au
            tpPh_jd.append(float(el['Tp_jd'])) # time of pericenter passage, Julian date
            WPh_deg.append(np.mod(float(el['Omega']),360)) # degrees
            wPh_deg.append(np.mod(float(el['w']),360)) # degrees
            iPh_deg.append(float(el['incl'])) # degrees
            # No clones for Pluto.
            dictionary = {'ePh':ePh,'qPh_au':qPh_au,'tpPh_jd':tpPh_jd,'WPh_deg':WPh_deg,'wPh_deg':wPh_deg,'iPh_deg':iPh_deg}
            df_out = pd.DataFrame(dictionary)
            df_out.to_csv(clones_file,index=False)
        else:
            # Now parse JSON files for covariance information.
            infile = 'ssd_json_' + des + '.json'
            with open(infile) as json_file:
                data = json.load(json_file)
            cov_array = data['orbit']['covariance']['data'] # orbels in order e,q,tp,node,peri,i
            # Row-major or column-major interpretation of covariance matrix doesn't matter,
            # because covariance matrix is symmetric.
            cov_array_2 = np.array([\
                [float(cov_array[0][0]),float(cov_array[0][1]),float(cov_array[0][2]),\
                     float(cov_array[0][3]),float(cov_array[0][4]),float(cov_array[0][5])],\
                [float(cov_array[1][0]),float(cov_array[1][1]),float(cov_array[1][2]),\
                     float(cov_array[1][3]),float(cov_array[1][4]),float(cov_array[1][5])],\
                [float(cov_array[2][0]),float(cov_array[2][1]),float(cov_array[2][2]),\
                     float(cov_array[2][3]),float(cov_array[2][4]),float(cov_array[2][5])],\
                [float(cov_array[3][0]),float(cov_array[3][1]),float(cov_array[3][2]),\
                     float(cov_array[3][3]),float(cov_array[3][4]),float(cov_array[3][5])],\
                [float(cov_array[4][0]),float(cov_array[4][1]),float(cov_array[4][2]),\
                     float(cov_array[4][3]),float(cov_array[4][4]),float(cov_array[4][5])],\
                [float(cov_array[5][0]),float(cov_array[5][1]),float(cov_array[5][2]),\
                     float(cov_array[5][3]),float(cov_array[5][4]),float(cov_array[5][5])] ])
            # Covariance and elements are specified at this epoch (jd).
            cov_epoch = data['orbit']['covariance']['epoch']
            cov_epoch = float(cov_epoch)
            # Save heliocentric orbital elements of the planets at the covariance epoch.
            epochP.append(cov_epoch)
            obj = Horizons(id='5',location=center,epochs=cov_epoch) # Jupiter barycenter
            el = obj.elements()
            eJh.append(float(el['e']))
            qJh_au.append(float(el['q']))
            tpJh_jd.append(float(el['Tp_jd']))
            WJh_deg.append(np.mod(float(el['Omega']),360))
            wJh_deg.append(np.mod(float(el['w']),360))
            iJh_deg.append(float(el['incl']))
            obj = Horizons(id='6',location=center,epochs=cov_epoch) # Saturn barycenter
            el = obj.elements()
            eSh.append(float(el['e']))
            qSh_au.append(float(el['q']))
            tpSh_jd.append(float(el['Tp_jd']))
            WSh_deg.append(np.mod(float(el['Omega']),360))
            wSh_deg.append(np.mod(float(el['w']),360))
            iSh_deg.append(float(el['incl']))
            obj = Horizons(id='7',location=center,epochs=cov_epoch) # Uranus barycenter
            el = obj.elements()
            eUh.append(float(el['e']))
            qUh_au.append(float(el['q']))
            tpUh_jd.append(float(el['Tp_jd']))
            WUh_deg.append(np.mod(float(el['Omega']),360))
            wUh_deg.append(np.mod(float(el['w']),360))
            iUh_deg.append(float(el['incl']))
            obj = Horizons(id='8',location=center,epochs=cov_epoch) # Neptune barycenter
            el = obj.elements()
            eNh.append(float(el['e']))
            qNh_au.append(float(el['q']))
            tpNh_jd.append(float(el['Tp_jd']))
            WNh_deg.append(np.mod(float(el['Omega']),360))
            wNh_deg.append(np.mod(float(el['w']),360))
            iNh_deg.append(float(el['incl']))
            dictionary = {'epochP_jd':epochP,
                  'eJh':eJh,'qJh_au':qJh_au,'tpJh_jd':tpJh_jd,'WJh_deg':WJh_deg,'wJh_deg':wJh_deg,'iJh_deg':iJh_deg,\
                  'eSh':eSh,'qSh_au':qSh_au,'tpSh_jd':tpSh_jd,'WSh_deg':WSh_deg,'wSh_deg':wSh_deg,'iSh_deg':iSh_deg,\
                  'eUh':eUh,'qUh_au':qUh_au,'tpUh_jd':tpUh_jd,'WUh_deg':WUh_deg,'wUh_deg':wUh_deg,'iUh_deg':iUh_deg,\
                  'eNh':eNh,'qNh_au':qNh_au,'tpNh_jd':tpNh_jd,'WNh_deg':WNh_deg,'wNh_deg':wNh_deg,'iNh_deg':iNh_deg}
            df_out = pd.DataFrame(dictionary)
            df_out.to_csv(planets_file,index=False)
            # Now retrieve the nominal orbital elements at epoch.
            e_nom  = data['orbit']['covariance']['elements'][0]['value']
            q_nom  = data['orbit']['covariance']['elements'][1]['value']
            tp_nom = data['orbit']['covariance']['elements'][2]['value']
            W_nom  = data['orbit']['covariance']['elements'][3]['value']
            w_nom  = data['orbit']['covariance']['elements'][4]['value']
            i_nom  = data['orbit']['covariance']['elements'][5]['value']
            e_nom  = float(e_nom)
            q_nom  = float(q_nom) # au
            tp_nom = float(tp_nom) # time of pericenter passage, Julian date
            W_nom  = float(W_nom) # degrees
            w_nom  = float(w_nom) # degrees
            i_nom  = float(i_nom) # degrees
            # The first line of clones_[des].csv is always the heliocentric orbels
            # of the Plutino and the planets at the epoch of the covariance matrix.
            ePh.append(e_nom)
            qPh_au.append(q_nom)
            tpPh_jd.append(tp_nom)
            WPh_deg.append(W_nom)
            wPh_deg.append(w_nom)
            iPh_deg.append(i_nom)
            epochP.append(cov_epoch)
            # Make nominal state vector before generating clones.
            state_nom = np.array([e_nom,q_nom,tp_nom,W_nom,w_nom,i_nom])
            # Generate clones.
            for j in range(Nclones):
                L = np.linalg.cholesky(cov_array_2)
                clone_1 = np.random.normal(0,1,6) # 6-element random Gaussian mean 0 stdev 1
                clone_2 = np.dot(L,clone_1) # Gaussian perturbations from zero correlated according to cov matrix
                clone_3 = clone_2 + state_nom # add perturbation to nominal orbit
                e_clone = clone_3[0]
                q_clone = clone_3[1]
                tp_clone = clone_3[2]
                W_clone = clone_3[3]
                w_clone = clone_3[4]
                i_clone = clone_3[5]
                if (e_clone<=0 or q_clone<=0): # don't want unphysical orbits
                    print('clone problems with ',iobj,des)
                W_clone = np.mod(W_clone,360)
                w_clone = np.mod(w_clone,360)
                # Record heliocentric orbels of Plutino.
                ePh.append(e_clone)
                qPh_au.append(q_clone)
                tpPh_jd.append(tp_clone)
                WPh_deg.append(W_clone)
                wPh_deg.append(w_clone)
                iPh_deg.append(i_clone)
        dictionary = {'ePh':ePh,'qPh_au':qPh_au,'tpPh_jd':tpPh_jd,'WPh_deg':WPh_deg,'wPh_deg':wPh_deg,'iPh_deg':iPh_deg}
        df_out = pd.DataFrame(dictionary)
        df_out.to_csv(clones_file,index=False)
    return
#%%
def get_GMdict():
# Convenient storage of outer planet masses.
# reference:
    GM_Sun = 1
    GM_Mercury = 1/6023600
    GM_Venus = 1/408523.71
    GM_EarthMoon = 1/328900.56
    GM_mars = 1/3098708
    GM_Jupiter = 1/1047.3486
    GM_Saturn = 1/3497.898
    GM_Uranus = 1/22902.98
    GM_Neptune = 1/19412.24
    # Add inner planet masses to Sun and renormalize outer planet masses.
    GM_Sun = GM_Sun + GM_Mercury + GM_Venus + GM_EarthMoon + GM_mars
    GM_Jupiter = GM_Jupiter / GM_Sun
    GM_Saturn  = GM_Saturn  / GM_Sun
    GM_Uranus  = GM_Uranus  / GM_Sun
    GM_Neptune = GM_Neptune / GM_Sun
    GM_Sun =     GM_Sun     / GM_Sun # Trivial, but in here for completeness.
    GMdict = {'Sun':GM_Sun,'Jupiter':GM_Jupiter,'Saturn':GM_Saturn,\
              'Uranus':GM_Uranus,'Neptune':GM_Neptune}
    return GMdict
#%%
def get_gmv08settings():
# Parameters for Kuiper belt objects that are worth classifying, according to
# 'Nomenclature in the Outer Solar System', Gladman/Marsden/VanLaerhoven 2008
# (gmv08), chapter 43 of 'The Solar System Beyond Neptune', U. of Arizona Press.
    tjmax = 3.05 # Jupiter Tisserand parameter used in gmv08 scheme.
    qmin = 7.35 # Minimum perihelion used in gmv08 scheme, au.
    a_Oort = 2000 # Oort cloud distance used in gmv08 scheme.
    e_Detached = 0.24 # Min eccentricity for Detached objects in gmv08 scheme.
    a_outer = 47.8 # Boundary for outer Classical belt in gmv08 scheme, au.
    a_inner = 39.4 # Boundary for inner Classical belt in gmv08 scheme, au.
    Nopp_min = 3 # Minimum number of oppositions to observe to pin down the orbit.
    gmv08settings = {'tjmax':tjmax,'qmin':qmin,'a_Oort':a_Oort,'e_Detached':e_Detached,\
                     'a_outer':a_outer,'a_inner':a_inner,'Nopp_min':Nopp_min}
    return gmv08settings
#%%
def get_sv20classifier():
# Machine learning classifier for Kuiper belt objects. Code is unchanged from
# code distributed with 'Machine learning classification of Kuiper belt populations'.
# Smullen/Volk 2020 (sv20), MNRAS 497:2, September 2020, pg 1391-1403,
# https://doi.org/10.1093/mnras/staa1935
# Code is found at https://github.com/rsmullen/KBO_Classifier
    import pandas as pd
    from sklearn.model_selection import train_test_split
    from sklearn.ensemble import GradientBoostingClassifier
    training_file = '00_KBO_features.csv'
    all_KBOs = pd.read_csv(training_file,skipinitialspace=True)
    secure_KBOs = all_KBOs[all_KBOs['Securely Classified']==True]
    all_types = list(set(secure_KBOs['Class']))
    types_dict = { all_types[i] : i for i in range( len(all_types) ) }
    int_dict = { i : all_types[i] for i in range( len(all_types) ) }
    classes = secure_KBOs['Class'].map(types_dict)
    features_train, features_test, classes_train, classes_test = train_test_split(secure_KBOs, classes, test_size=0.3, random_state=30)
    features_train.drop(['MPC ID', 'Securely Classified', 'Class'], axis=1, inplace=True)
    features_train = features_train.to_numpy()
    features_test.drop(['MPC ID', 'Securely Classified', 'Class'], axis=1, inplace=True)
    features_test = features_test.to_numpy()
    classifier = GradientBoostingClassifier( learning_rate=0.1, loss='deviance', max_depth=3, max_features='log2', n_estimators=130, random_state=30 )
    classifier.fit(features_train, classes_train)
    return classifier, int_dict, types_dict
#%% Integrate a single object and return orbels_dict for the Yu2018 resonance classification algorithm.
def integrate_single_tno(sim,tmax_yrs,tstep_yrs):
    import numpy as np
    sim.N_active = sim.N - 1 # tno is massless
    sim.move_to_com()
    tmax = tmax_yrs * 2*np.pi
    tstep = tstep_yrs * 2*np.pi
    t_list = []
    a_list = []
    a_Jupiter_list = []
    a_Neptune_list = []
    tj_list = []
    q_list = []
    e_list = []
    inc_list = []
    w_list = []
    W_list = []
    M_list = []
    l_list = []
    l_Neptune_list = []
    pomega_list = []
    pomega_Neptune_list = []
    f_list = []
    f_Neptune_list = []
    P_list = []
    P_Neptune_list = []
    theta_list = []
    theta_Neptune_list = []
    sim.t = 0
    t = 0
    com = sim.calculate_com()
    t_list.append(t)
    Jupiter = sim.particles['Jupiter']
    Neptune = sim.particles['Neptune']
    tno = sim.particles['Plutino']
    orbit_Jupiter = Jupiter.calculate_orbit(primary=com)
    orbit_Neptune = Neptune.calculate_orbit(primary=com)
    orbit_tno = tno.calculate_orbit(primary=com)
    a = orbit_tno.a
    a_Jupiter = orbit_Jupiter.a
    a_Neptune = orbit_Neptune.a
    e = orbit_tno.e
    q = a * (1-e)
    inc = orbit_tno.inc
    tj = Tisserand(a,a_Jupiter,inc,e)
    w = orbit_tno.omega
    W = orbit_tno.Omega
    M = orbit_tno.M
    l = orbit_tno.l
    l_Neptune = orbit_Neptune.l
    pomega = orbit_tno.pomega
    pomega_Neptune = orbit_Neptune.pomega
    f = orbit_tno.f
    f_Neptune = orbit_Neptune.f
    P = orbit_tno.P
    P_Neptune = orbit_Neptune.P
    theta = orbit_tno.theta
    theta_Neptune = orbit_Neptune.theta
    a_list.append(a)
    a_Jupiter_list.append(a_Jupiter)
    a_Neptune_list.append(a_Neptune)
    q_list.append(q)
    e_list.append(e)
    inc_list.append(inc)
    tj_list.append(tj)
    w_list.append(w)
    W_list.append(W)
    M_list.append(M)
    l_list.append(l)
    l_Neptune_list.append(l_Neptune)
    pomega_list.append(pomega)
    pomega_Neptune_list.append(pomega_Neptune)
    f_list.append(f)
    f_Neptune_list.append(f_Neptune)
    P_list.append(P)
    P_Neptune_list.append(P_Neptune)
    theta_list.append(theta)
    theta_Neptune_list.append(theta_Neptune)
    while t < tmax:
        t = t + tstep
        sim.integrate(t,exact_finish_time=1)
        # # record stuff
        com = sim.calculate_com()
        t_list.append(t)
        Jupiter = sim.particles['Jupiter']
        Neptune = sim.particles['Neptune']
        tno = sim.particles['Plutino'] # 0sun,1jup,2sat,3ur,4nep
        orbit_Jupiter = Jupiter.calculate_orbit(primary=com)
        orbit_Neptune = Neptune.calculate_orbit(primary=com)
        orbit_tno = tno.calculate_orbit(primary=com)
        a = orbit_tno.a
        a_Jupiter = orbit_Jupiter.a
        a_Neptune = orbit_Neptune.a
        e = orbit_tno.e
        q = a * (1-e)
        inc = orbit_tno.inc
        tj = Tisserand(a,a_Jupiter,inc,e)
        w = orbit_tno.omega
        W = orbit_tno.Omega
        M = orbit_tno.M
        l = orbit_tno.l
        l_Neptune = orbit_Neptune.l
        pomega = orbit_tno.pomega
        pomega_Neptune = orbit_Neptune.pomega
        f = orbit_tno.f
        f_Neptune = orbit_Neptune.f
        P = orbit_tno.P
        P_Neptune = orbit_Neptune.P
        theta = orbit_tno.theta
        theta_Neptune = orbit_Neptune.theta
        a_list.append(a)
        a_Jupiter_list.append(a_Jupiter)
        a_Neptune_list.append(a_Neptune)
        q_list.append(q)
        e_list.append(e)
        inc_list.append(inc)
        tj_list.append(tj)
        w_list.append(w)
        W_list.append(W)
        M_list.append(M)
        l_list.append(l)
        l_Neptune_list.append(l_Neptune)
        pomega_list.append(pomega)
        pomega_Neptune_list.append(pomega_Neptune)
        f_list.append(f)
        f_Neptune_list.append(f_Neptune)
        P_list.append(P)
        P_Neptune_list.append(P_Neptune)
        theta_list.append(theta)
        theta_Neptune_list.append(theta_Neptune)
    # # done integrating, now save stuff to a dictionary
    t_array = np.array(t_list)
    t_yrs_array = t_array/2/np.pi
    a_array = np.array(a_list)
    a_Jupiter_array = np.array(a_Jupiter_list)
    a_Neptune_array = np.array(a_Neptune_list)
    tj_array = np.array(tj_list)
    q_array = np.array(q_list)
    e_array = np.array(e_list)
    inc_array = np.array(inc_list)
    w_array = np.array(w_list)
    W_array = np.array(W_list)
    M_array = np.array(M_list)
    l_array = np.array(l_list)
    l_Neptune_array = np.array(l_Neptune_list)
    pomega_array = np.array(pomega_list)
    pomega_Neptune_array = np.array(pomega_Neptune_list)
    f_array = np.array(f_list)
    f_Neptune_array = np.array(f_Neptune_list)
    P_array = np.array(P_list)
    P_Neptune_array = np.array(P_Neptune_list)
    theta_array = np.array(theta_list)
    theta_Neptune_array = np.array(theta_Neptune_list)
    w_array = np.mod(w_array,2*np.pi)
    W_array = np.mod(W_array,2*np.pi)
    M_array = np.mod(M_array,2*np.pi)
    l_array = np.mod(l_array,2*np.pi)
    l_Neptune_array = np.mod(l_Neptune_array,2*np.pi)
    pomega_array = np.mod(pomega_array,2*np.pi)
    pomega_Neptune_array = np.mod(pomega_Neptune_array,2*np.pi)
    f_array = np.mod(f_array,2*np.pi)
    f_Neptune_array = np.mod(f_Neptune_array,2*np.pi)
    theta_array = np.mod(theta_array,2*np.pi)
    theta_Neptune_array = np.mod(theta_Neptune_array,2*np.pi)
    inc_array = np.degrees(inc_array)
    w_array = np.degrees(w_array)
    W_array = np.degrees(W_array)
    M_array = np.degrees(M_array)
    l_array = np.degrees(l_array)
    l_Neptune_array = np.degrees(l_Neptune_array)
    pomega_array = np.degrees(pomega_array)
    pomega_Neptune_array = np.degrees(pomega_Neptune_array)
    f_array = np.degrees(f_array)
    f_Neptune_array = np.degrees(f_Neptune_array)
    P_yrs_array = P_array/2/np.pi
    P_yrs_Neptune_array = P_Neptune_array/2/np.pi
    theta_array = np.degrees(theta_array)
    theta_Neptune_array = np.degrees(theta_Neptune_array)
    orbels_dict = {'t':t_yrs_array,'a':a_array,'a_Jupiter':a_Jupiter_array,\
        'a_Neptune':a_Neptune_array,'tj':tj_array,'q':q_array,'e':e_array,\
        'inc':inc_array,'w':w_array,'W':W_array,'M':M_array,'l':l_array,\
        'l_Neptune':l_Neptune_array,'pomega':pomega_array,'pomega_Neptune':\
        pomega_Neptune_array,'f':f_array,'f_Neptune':f_Neptune_array,\
        'P':P_yrs_array,'P_Neptune':P_yrs_Neptune_array,'theta':theta_array,\
        'theta_Neptune':theta_Neptune_array}
    sim = None
    return orbels_dict
#%%
def inverse_transform_sampling(data, n_bins=40, n_samples=1000):
    # also from Stack Exchange or a Google search, I forgot to copy the link
    import scipy.interpolate as interpolate
    import numpy as np
    hist, bin_edges = np.histogram(data, bins=n_bins, density=True)
    cum_values = np.zeros(bin_edges.shape)
    cum_values[1:] = np.cumsum(hist*np.diff(bin_edges))
    inv_cdf = interpolate.interp1d(cum_values, bin_edges)
    r = np.random.rand(n_samples)
    return inv_cdf(r)
#%%
def JDfun(year,month,day,hour,minute,second):
# Compute Julian date from calendar date.
# From Vallado 'Fundamentals of Astrodynamics and Applications' 4th ed. pg 183
# Valid for years 1900 to 2100, CE.
    import numpy as np
    JD = 367*year - np.floor(7/4*(year+np.floor(1/12*(month+9)))) + \
        np.floor(275*month/9) + day + 1721013.5 + 1/24*(hour+1/60*(second/60+minute))
    return JD
#%% make list of mmrs to check
def make_mmrs(amin,amax,a_neptune,bigmax,smallmax,ordermax):
    import math
    import pandas as pd
    big_list = []
    small_list = []
    period_ratio_list = []
    sma_list = []
    for big in range(bigmax):
        for small in range(smallmax):
            big = big + 1
            small = small + 1 # adjust for Python's 0-indexing
            gcd = math.gcd(big,small)
            atno = (big/small)**(2/3) * a_neptune
            if gcd == 1 and (amin <= atno < amax) and (big>=small) and (big-small<=ordermax):
                big_list.append(big)
                small_list.append(small)
                period_ratio_list.append(big/small)
                sma_list.append(atno)
    data = {'big':big_list,\
            'small':small_list,\
            'period_ratio':period_ratio_list,\
            'semimajor_axis':sma_list}
    df_mmrs = pd.DataFrame.from_dict(data)
    return df_mmrs
#%%
def parsedata_sv20classifier(data,classifier,int_dict):
# This function computes the necessary features to classify a KBO.
# Data MUST be a 101 row x 6 column array.
# Columns are t, a, e, i, Omega, omega.
# Rows are different time outputs: MUST be 1000yr outputs, ie [0, 1E3, 2E3....99E3,100E3].
# Returns features for classification.
# This code is unchanged from sv20, ie
# Smullen, Rachel A., and Kathryn Volk.
# 'Machine learning classification of Kuiper belt populations.'
# Monthly Notices of the Royal Astronomical Society 497.2 (2020): 1391-1403.
# The classification simulation and data parsing is copied from sv20 code with
# minimal changes to not have to query Horizons through Rebound when setting up a sim,
# because Rebound's Horizons query is very slow.
    import numpy as np
    # Take stats of simulations.
    initials = data[0,1:] # a, e, i, Omega, omega
    finals = data[-1,1:]
    mins = np.amin(data[:,1:],axis = 0)
    maxes = np.amax(data[:,1:],axis = 0)
    dels = maxes-mins
    means = np.mean(data[:,1:],axis = 0)
    stdev = np.std(data[:,1:],axis = 0)
    # Take time derivatives.
    diffs = data[1:,:]-data[:-1,:]
    dxdt = diffs[:,1:]/diffs[:,0, np.newaxis] # Add on new axis to time to give same dimensionality as the numerator.
    mindxdt = np.amin(dxdt,axis = 0)
    meandxdt = np.mean(dxdt,axis = 0)
    maxdxdt = np.amax(dxdt,axis = 0)
    deldxdt = maxdxdt-mindxdt
    # Rearrange data into the order we want.
    arrs = [initials,finals,mins,means,maxes,stdev,dels,mindxdt,meandxdt,maxdxdt,deldxdt]
    inds = [0,1,2,3,4] # a, e, i, Omega, omega
    features = []
    ## Features contains all x values, then all y, etc: xi, xf, xmin, xmean, xmax, xsigma, Deltax, xdotmin, xdotmean, xdotmax.
    for i in inds:
        for a in arrs:
            features += [a[i]]
    features_out = np.array(features).reshape(1,-1) # Make sure features is a 2d array.
    prediction = classifier.predict_proba(features_out) # Predict the probabilities of class membership for object.
    if np.max(prediction) == prediction[0][0]:
        category = int_dict[0]
    elif np.max(prediction) == prediction[0][1]:
        category = int_dict[1]
    elif np.max(prediction) == prediction[0][2]:
        category = int_dict[2]
    elif np.max(prediction) == prediction[0][3]:
        category = int_dict[3]
    return category,prediction
#%%
def rand_t_marginal(kappa,p,N=1):
    """
    https://dlwhittenbury.github.io/ds-2-sampling-and-visualising-the-von-mises-fisher-distribution-in-p-dimensions.html
        rand_t_marginal(kappa,p,N=1)
        ============================
        Samples the marginal distribution of t using rejection sampling of Wood [3].
        INPUT:
            * kappa (float) - concentration
            * p (int) - The dimension of the generated samples on the (p-1)-dimensional hypersphere.
                - p = 2 for the unit circle $\mathbb{S}^{1}$
                - p = 3 for the unit sphere $\mathbb{S}^{2}$
            Note that the (p-1)-dimensional hypersphere $\mathbb{S}^{p-1} \subset \mathbb{R}^{p}$ and the
            samples are unit vectors in $\mathbb{R}^{p}$ that lie on the sphere $\mathbb{S}^{p-1}$.
            * N (int) - number of samples
        OUTPUT:
            * samples (array of floats of shape (N,1)) - samples of the marginal distribution of t
    """
    import numpy as np
    # Check kappa >= 0 is numeric
    # if (kappa < 0) or ((type(kappa) is not float) and (type(kappa) is not int)):
    if (kappa < 0):
        raise Exception("kappa must be a non-negative number.")
    if (p<=0) or (type(p) is not int):
        raise Exception("p must be a positive integer.")
    # Check N>0 and is an int
    if (N<=0) or (type(N) is not int):
        raise Exception("N must be a non-zero positive integer.")
    # Start of algorithm
    b = (p - 1.0) / (2.0 * kappa + np.sqrt(4.0 * kappa**2 + (p - 1.0)**2 ))
    x0 = (1.0 - b) / (1.0 + b)
    c = kappa * x0 + (p - 1.0) * np.log(1.0 - x0**2)
    samples = np.zeros((N,1))
    # Loop over number of samples
    for i in range(N):
        # Continue unil you have an acceptable sample
        while True:
            # Sample Beta distribution
            Z = np.random.beta( (p - 1.0)/2.0, (p - 1.0)/2.0 )
            # Sample Uniform distribution
            U = np.random.uniform(low=0.0,high=1.0)
            # W is essentially t
            W = (1.0 - (1.0 + b) * Z) / (1.0 - (1.0 - b) * Z)
            # Check whether to accept or reject
            if kappa * W + (p - 1.0)*np.log(1.0 - x0*W) - c >= np.log(U):
                # Accept sample
                samples[i] = W
                break
    return samples
#%%
def rand_uniform_hypersphere(N,p):
    """
    https://dlwhittenbury.github.io/ds-2-sampling-and-visualising-the-von-mises-fisher-distribution-in-p-dimensions.html
        rand_uniform_hypersphere(N,p)
        =============================
        Generate random samples from the uniform distribution on the (p-1)-dimensional
        hypersphere $\mathbb{S}^{p-1} \subset \mathbb{R}^{p}$. We use the method by
        Muller [1], see also Ref. [2] for other methods.
        INPUT:
            * N (int) - Number of samples
            * p (int) - The dimension of the generated samples on the (p-1)-dimensional hypersphere.
                - p = 2 for the unit circle $\mathbb{S}^{1}$
                - p = 3 for the unit sphere $\mathbb{S}^{2}$
            Note that the (p-1)-dimensional hypersphere $\mathbb{S}^{p-1} \subset \mathbb{R}^{p}$ and the
            samples are unit vectors in $\mathbb{R}^{p}$ that lie on the sphere $\mathbb{S}^{p-1}$.
    References:
    [1] Muller, M. E. "A Note on a Method for Generating Points Uniformly on N-Dimensional Spheres."
    Comm. Assoc. Comput. Mach. 2, 19-20, Apr. 1959.
    [2] https://mathworld.wolfram.com/SpherePointPicking.html
    """
    import numpy as np
    if (p<=0) or (type(p) is not int):
        raise Exception("p must be a positive integer.")
    # Check N>0 and is an int
    if (N<=0) or (type(N) is not int):
        raise Exception("N must be a non-zero positive integer.")
    v = np.random.normal(0,1,(N,p))
    v = np.divide(v,np.linalg.norm(v,axis=1,keepdims=True))
    return v
#%%
def rand_von_mises_fisher(mu,kappa,N=1):
    """
    https://dlwhittenbury.github.io/ds-2-sampling-and-visualising-the-von-mises-fisher-distribution-in-p-dimensions.html
        rand_von_mises_fisher(mu,kappa,N=1)
        ===================================
        Samples the von Mises-Fisher distribution with mean direction mu and concentration kappa.
        INPUT:
            * mu (array of floats of shape (p,1)) - mean direction. This should be a unit vector.
            * kappa (float) - concentration.
            * N (int) - Number of samples.
        OUTPUT:
            * samples (array of floats of shape (N,p)) - samples of the von Mises-Fisher distribution
            with mean direction mu and concentration kappa.
    """
    import numpy as np
    import numpy.matlib
    from scipy.linalg import null_space
    # Check that mu is a unit vector
    eps = 10**(-8) # Precision
    norm_mu = np.linalg.norm(mu)
    if abs(norm_mu - 1.0) > eps:
        raise Exception("mu must be a unit vector.")
    # Check kappa >= 0 is numeric
    # if (kappa < 0) or ((type(kappa) is not float) and (type(kappa) is not int)):
    if (kappa < 0):
        raise Exception("kappa must be a non-negative number.")
    # Check N>0 and is an int
    if (N<=0) or (type(N) is not int):
        raise Exception("N must be a non-zero positive integer.")
    # Dimension p
    p = len(mu)
    # Make sure that mu has a shape of px1
    mu = np.reshape(mu,(p,1))
    # Array to store samples
    samples = np.zeros((N,p))
    #  Component in the direction of mu (Nx1)
    t = rand_t_marginal(kappa,p,N)
    # Component orthogonal to mu (Nx(p-1))
    xi = rand_uniform_hypersphere(N,p-1)
    # von-Mises-Fisher samples Nxp
    # Component in the direction of mu (Nx1).
    # Note that here we are choosing an
    # intermediate mu = [1, 0, 0, 0, ..., 0] later
    # we rotate to the desired mu below
    samples[:,[0]] = t
    # Component orthogonal to mu (Nx(p-1))
    samples[:,1:] = np.matlib.repmat(np.sqrt(1 - t**2), 1, p-1) * xi
    # Rotation of samples to desired mu
    O = null_space(mu.T)
    R = np.concatenate((mu,O),axis=1)
    samples = np.dot(R,samples.T).T
    return samples
#%%
def rodrigues_rotation(vold,k,theta):
# Rodrigues' rotation formula
# vold is a vector in R3
# k is a unit vector describing an axis of rotation about which vold rotates
# theta is the angle of rotation according to the right hand rule
    import numpy as np
    # equation A2 in the paper
    term1 = vold * np.cos(theta)
    term2 = np.cross(k,vold) * np.sin(theta)
    term3 = k * np.dot(k,vold) * (1-np.cos(theta))
    vnew = term1 + term2 + term3
    return vnew
#%%
def sbdb_reduce(JD,datestr):
# Reduce MPC and JPL SBDB files to a more manageable number of candidate Plutinos.
    import pandas as pd
    import numpy as np
    from astroquery.jplhorizons import Horizons
    # Get gmv08 criteria for which objects we can try to classify.
    gmv08settings = get_gmv08settings()
    tjmax = gmv08settings['tjmax']
    qmin = gmv08settings['qmin']
    max_fractional_sigma_a = 0.05 # This rule of thumb is from vm17,
    # Volk/Malhotra 2017, 'The curiously warped mean plane of the Kuiper belt'.
    mpcfile = '00_MPCORB_reduced_' + datestr + '.csv'
    jplfile = '00_sbdb_query_results_' + datestr + '.csv'
    outfile = 'sbdb_reduce_output.csv'
    # Read in objects in the reduced MPC database.
    dfmpc = pd.read_csv(mpcfile,low_memory=False)
    # Read in objects in the SBDB database.
    dfjpl = pd.read_csv(jplfile,low_memory=False)
    Nstart = dfjpl.shape[0]
    # Drop objects in the JPL list with too much semimajor axis uncertainty to be worth classifying.
    heliocentric_a_list = dfjpl['a'].tolist()
    temp_sigma_a_list = dfjpl['sigma_a'].tolist()
    for iobj in range(Nstart):
        if np.isnan(temp_sigma_a_list[iobj]):
            temp_sigma_a_list[iobj] = 0 # Need this to make sure we keep Pluto.
    dfjpl['sigma_a'] = temp_sigma_a_list
    heliocentric_sigma_a_list = dfjpl['sigma_a'].tolist()
    heliocentric_a_array = np.array(heliocentric_a_list)
    heliocentric_sigma_a_array = np.array(heliocentric_sigma_a_list)
    heliocentric_fractional_sigma_a_array = heliocentric_sigma_a_array / heliocentric_a_array
    dfjpl['fractional_sigma_a'] = heliocentric_fractional_sigma_a_array
    dfjpl = dfjpl[dfjpl['fractional_sigma_a']<=max_fractional_sigma_a]
    # Now we find objects that are in BOTH the trimmed MPC list and the trimmed JPL list.
    # This ensures that opposition count is high enough and fractional semimajor axis
    # uncertainty is low enough.
    mpcdes = dfmpc['packed_designation'].tolist()
    jpldes = dfjpl['pdes'].tolist()
    full_name = dfjpl['full_name'].tolist()
    commondes = []
    Njpl = len(jpldes)
    Nmpc = len(mpcdes)
    for i in range(Njpl):
        des = jpldes[i]
        fullname = full_name[i]
        des = des.lstrip()
        des = des.rstrip()
        fullname = fullname.lstrip()
        fullname = fullname.rstrip()
        if fullname.startswith('C/'):
            2-2 # Comet; don't add it to the list of objects in common; we don't want it.
        else:
            des2 = sbdb_des_to_packed_mpc_des(des)
            if des2 in mpcdes:
                commondes.append(des2)
    for i in range(Nmpc):
        des = mpcdes[i]
        if (des in jpldes) and (des not in commondes):
            commondes.append(des)
    Ncom = len(commondes)
    # Retrieve semimajor axis of Neptune.
    center = '500@0' # Solar System barycenter
    obj = Horizons(id='8',location=center,epochs=JD)
    el = obj.elements()
    a_Neptune = float(el['a']) # au
    # Now we retrieve barycentric elements for each object.
    a_list = []
    e_list = []
    i_list = []
    for i in range(Ncom):
        des = commondes[i]
        unpacked = unpack(des)
        if des == 'D4340':
            unpacked = '9' # We want the Pluto-Charon barycenter, not the Pluto body center.
        obj = Horizons(id=unpacked,location=center,epochs=JD)
        el = obj.elements()
        a = float(el['a']) # au
        e = float(el['e'])
        incl = float(el['incl']) # deg
        incl = np.mod(incl,360)
        a_list.append(a)
        e_list.append(e)
        i_list.append(incl)
    orbels_dict = {'packed_designation':commondes,'a_au':a_list,'e':e_list,'i_deg':i_list}
    dfcom = pd.DataFrame.from_dict(orbels_dict)
    # Now we do a fine cut to barycentric semimajor axis of 34 to 150 au.
    # We also do a perihelion cut and Tisserand cut to match gmv08 resonance eligibility requirements.
    a_array = np.array(a_list)
    e_array = np.array(e_list)
    q_array = a_array * (1-e_array)
    name = '5' # Jupiter
    obj = Horizons(id=name,location=center,epochs=JD)
    el = obj.elements()
    a = float(el['a'])
    a_perturber = a
    elliptical_comet_list = []
    elliptical_comet_count = 0
    for i in range(Ncom):
        a = a_array[i]
        e = e_array[i]
        inc = np.radians(i_list[i])
        tj = Tisserand(a,a_perturber,inc,e)
        q = q_array[i]
        if q<qmin and tj<tjmax:
            elliptical_comet_list.append(2)
            elliptical_comet_count = elliptical_comet_count + 1
        else:
            elliptical_comet_list.append(0)
    dfcom['elliptical_comet'] = elliptical_comet_list
    dfcom = dfcom[dfcom['elliptical_comet']<1] # Drop elliptical comets.
    dfcom = dfcom[dfcom['e']<1] # Drop hyperbolic comets.
    # Semimajor axis limits are still a rough cut to make sure we don't miss any Plutinos.
    dfcom = dfcom[dfcom['a_au']>=1.2*a_Neptune] # Enforce semimajor axis limits.
    dfcom = dfcom[dfcom['a_au']<=1.4*a_Neptune] # Enforce semimajor axis limits.
    dfcom = dfcom.drop(columns=['elliptical_comet']) # We don't care about saving this info.
    dfcom = dfcom.drop(columns=['a_au']) # We don't care about saving this info.
    dfcom = dfcom.drop(columns=['e']) # We don't care about saving this info.
    dfcom = dfcom.drop(columns=['i_deg']) # We don't care about saving this info.
    Ncom = dfcom.shape[0]
    dfcom.to_csv(outfile,index=False)
    return Nstart,Ncom
#%%
def sbdb_des_to_packed_mpc_des(des1):
# Bring in SBDB designations, get out MPC designations.
# See Minor Planet Center information on packed designations.
    from bidict import bidict
    onedigit_bidict = bidict({'0':'0','1':'1','2':'2','3':'3','4':'4','5':'5','6':'6',\
        '7':'7','8':'8','9':'9','A':'10','B':'11','C':'12','D':'13','E':'14',\
        'F':'15','G':'16','H':'17','I':'18','J':'19','K':'20','L':'21','M':'22',\
        'N':'23','O':'24','P':'25','Q':'26','R':'27','S':'28','T':'29','U':'30',\
        'V':'31','W':'32','X':'33','Y':'34','Z':'35','a':'36','b':'37','c':'38',\
        'd':'39','e':'40','f':'41','g':'42','h':'43','i':'44','j':'45','k':'46',\
        'l':'47','m':'48','n':'49','o':'50','p':'51','q':'52','r':'53','s':'54',\
        't':'55','u':'56','v':'57','w':'58','x':'59','y':'60','z':'61'})
    century_bidict = bidict({'I':'18','J':'19','K':'20'})
    des1 = str(des1)
    if len(des1) == 5:
        des2 = des1 # leave it unchanged, example is 15760 for Albion
    if len(des1) == 6: # change the first two numerical digits to a letter
        first_two_digits = des1[0:2]
        last_four_digits = des1[2:6]
        one_digit = onedigit_bidict.inverse[first_two_digits]
        des2 = one_digit + last_four_digits
    if len(des1) > 6:
        century = des1[0:2]
        year = des1[2:4]
        space = des1[4]
        letter1 = des1[5]
        letter2 = des1[6]
        if len(des1) == 7:
            number = 0
        elif len(des1) == 8:
            number = des1[7]
        elif len(des1) == 9:
            number = des1[7:9]
        elif len(des1) == 10:
            number = des1[7:10]
        else:
            raise Exception('Unrecognized sbdb designation format.',des1)
        number = int(number)
        if space != ' ':
            raise Exception('Unrecognized sbdb designation format.',des1)
        century = century_bidict.inverse[century]
        if number == 0:
            code = '00'
        elif number < 10:
            code = '0' + str(number)
        elif number < 100:
            code = str(number)
        elif number < 620: # 619 is the max number with the MPC packing scheme
            numberstr = str(number)
            first_two_digits = numberstr[0:2]
            last_digit = numberstr[2]
            one_digit = onedigit_bidict.inverse[first_two_digits]
            code = one_digit + last_digit
        else:
            raise Exception('Unrecognized sbdb designation format.',des1)
        des2 = century + year + letter1 + code + letter2
    return des2
#%%
def Tisserand(a,a_perturber,inc,ecc):
# Calculate Tisserand parameter with respect to a perturber. Inclination in radians.
    import numpy as np
    tp = a_perturber/a + 2 * np.cos(inc) * np.sqrt(a/a_perturber*(1-ecc*ecc))
    return tp
#%%
# def truncated_rayleigh_pdf(x, sigma):
def truncated_rayleigh_pdf(x,kappa):
    import numpy as np
    # equation 26 in the paper
    sigma = 1/np.sqrt(kappa)
    # equation 28 in the paper
    C_R = 1/(1-np.exp(-np.pi**2/2/sigma**2))
    # equation 27 in the paper
    return C_R/sigma**2 * x * np.exp(-x**2/2/sigma**2)
#%%
def unpack(designation):
# Unpack Minor Planet Center MPC ID.
    packed_designation = designation.lstrip()
    packed_designation = packed_designation.rstrip()
    Nt = len(packed_designation)
    onedigit_dict = {'0':'0','1':'1','2':'2','3':'3','4':'4','5':'5','6':'6',\
        '7':'7','8':'8','9':'9','A':'10','B':'11','C':'12','D':'13','E':'14',\
        'F':'15','G':'16','H':'17','I':'18','J':'19','K':'20','L':'21','M':'22',\
        'N':'23','O':'24','P':'25','Q':'26','R':'27','S':'28','T':'29','U':'30',\
        'V':'31','W':'32','X':'33','Y':'34','Z':'35','a':'36','b':'37','c':'38',\
        'd':'39','e':'40','f':'41','g':'42','h':'43','i':'44','j':'45','k':'46',\
        'l':'47','m':'48','n':'49','o':'50','p':'51','q':'52','r':'53','s':'54',\
        't':'55','u':'56','v':'57','w':'58','x':'59','y':'60','z':'61'}
    century_dict = {'I':'18','J':'19','K':'20'}
    if Nt == 4:
        unpacked_designation = packed_designation
    if Nt == 5:
        onedigit = packed_designation[0]
        fourdigits = packed_designation[1:5]
        if onedigit in onedigit_dict:
            onedigit = onedigit_dict[onedigit]
            unpacked_designation = onedigit + fourdigits
        else:
            unpacked_designation = '-99999'
    elif Nt == 7:
        century = packed_designation[0]
        year = packed_designation[1:3]
        letter1 = packed_designation[3]
        code1 = packed_designation[4]
        code2 = packed_designation[5]
        letter2 = packed_designation[6]
        space = ' '
        if century in century_dict:
            century = century_dict[century]
            if code1 == '0' and code2 == '0':
                code1 = ''
                code2 = ''
            else:
                if code1 in onedigit_dict:
                    code1 = onedigit_dict[code1]
                    if code1 == '0':
                        code1 = ''
            unpacked_designation = century + year + space + letter1 + letter2 + \
                code1 + code2
        else:
            unpacked_designation = '-99999'
    unpacked_designation = unpacked_designation.lstrip()
    unpacked_designation = unpacked_designation.rstrip()
    return unpacked_designation
