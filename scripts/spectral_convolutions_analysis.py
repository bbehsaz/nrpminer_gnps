# def generate_convolutions(
# 	standardAminoMasses,peaks,e,thresh,charge,representative, precursorMass):
# 	constituentMonoMers = []
# 	finalClusters, ALLcleanedAutoconvolutions, precursorMass, charge = generateAllAutconv(peaks,e,thresh,representative,standardAminoMasses,precursorMass,charge)
# 	return finalClusters, ALLcleanedAutoconvolutions, precursorMass, charge, standardAminoMasses, e

def median(lst):
    quotient, remainder = divmod(len(lst), 2)
    if remainder:
        return sorted(lst)[quotient]
    return sum(sorted(lst)[quotient - 1:quotient + 1]) / 2.


def generate_aa_convolutions_vectorized(
        standardAminoMasses, peaks, e, realPepMass,num_aa_threshold,peaks_masses,convolthresh):
    sattelites_removed_linked_cyclopep_clusters = {}
    num_aa = 0
    peptide = False
    convolthresh = 3
    for offset in sorted(standardAminoMasses):
        convol_num = convolution(peaks, peaks, offset, e, realPepMass,peaks_masses,convolthresh)
        if convol_num>convolthresh:
            num_aa += 1
            if num_aa == num_aa_threshold:
                peptide = True
                break
    return peptide


    # return sattelites_removed_linked_cyclopep_clusters
    # finalconvols = {}
    # clustermass_sattelites_removed_linked_cyclopep_clusters = {}
    # for aa in sattelites_removed_linked_cyclopep_clusters:
    #     sumofalldists = 0
    #     totalnumalldists = 0
    #     alldistances = []
    #     for group in sattelites_removed_linked_cyclopep_clusters[aa]:
    #         if len(sattelites_removed_linked_cyclopep_clusters[aa]) > 0:
    #             distances = [
    #                 min([abs(pair[1] - pair[0]) for pair in sattelites_removed_linked_cyclopep_clusters[aa][group]],
    #                     key=lambda x: abs(aa - x))]
    #             # distances = [abs(pair[1]-pair[0]) for pair in sattelites_removed_linked_cyclopep_clusters[aa][group]]
    #             alldistances += distances
    #             sumofalldists += sum(distances)
    #             totalnumalldists += len(distances)
    #     if totalnumalldists != 0:
    #         # clustermass = sumofalldists/(len(alldistances)*1.0)
    #         clustermass = median(alldistances)
    #     else:
    #         clustermass = 0
    #     if abs(clustermass - aa) < e:
    #         clustermass_sattelites_removed_linked_cyclopep_clusters[aa] = sattelites_removed_linked_cyclopep_clusters[
    #             aa]
    #     else:
    #         clustermass_sattelites_removed_linked_cyclopep_clusters[aa] = {}

    # return clustermass_sattelites_removed_linked_cyclopep_clusters


def convolution(spectrum1, spectrum2, offset, e, realPepMass,peaks,convolthresh):
    convolutionPairPeaks_list = []
    convolution_number = 0

    def make_binary_addE(spectrumVector):
        binary_spectrumVector_e = [0] * 2000000
        for p in range(len(spectrumVector)):
            if spectrumVector[p] > 0:
                for x in range(int(p - (5 * round(e, 3) * 1000)), int(p + (5 * round(e, 3) * 1000))):
                    binary_spectrumVector_e[x] = 1
        return binary_spectrumVector_e

    radius = int(5*e*1000)
    metrequirement = False #Keep track if this offset has already met the minimum requirements
    # for peak1 in range(min(int(realPepMass * 1000), len(spectrum1))-int(offset*1000)-1):
    for peak in peaks:
        if metrequirement:
            break
        peak1 = int(peak * 1000)
        peak2_center = int(peak1 + (offset * 1000))
        min_range=max(peak2_center - radius, 0)
        print "----"
        print "----"
        print min_range
        max_range=min(2000000, peak2_center + radius)
        print max_range
        for peak2 in range(min_range,max_range):
            print peak2
            if spectrum2[peak2]>0:
                convolution_number += 1
                if convolution_number>convolthresh:
                    metrequirement = True
                break  

        # if sum([1 for peak2 in range(min_range,max_range+1) if spectrum2[peak2]>0])>0:
        #     convolution_number += 1

        # print [1 for peak2 in range(min_range,max_range+1) if spectrum2[peak2]>0]
        # print len(range(min_range,max_range))
        # print convolution_number
        # for peak2 in range(max(peak2_center - radius, 0),
        #                    min(2000000, peak2_center-radius):
        # for peak2 in range(max(int(peak1 + (offset - 5 * e) * 1000), 0),
        #                    min(2000000, int(peak1 + (offset + 5 * e) * 1000) + 1)):
        #     if spectrum2[peak2] > 0:
        #         convolution_number += 1

                # convolutionPairPeaks_list.append(
                #     (round(float(peak1) / 1000.0, 3), round(float(peak2) / 1000.0, 3)))  # peakpair
                # break
    # print convolution_number
    # exit()

    return convolution_number
    # for peak1 in [peak for peak in range(min(int(realPepMass*1000), len(spectrum1))) if spectrum1[peak]>0]:
    # 	peak2 = peak1+int(round(offset,3)*1000)
    # 	if spectrum2_vector_5e[peak2]> 0:
    # 		convolutionPairPeaks_list.append(( round(float(peak1)/1000.0,3),round(float(peak2)/1000.0,3) ) ) #peakpair
    '''
    redundancy_removed = {}

    current_group = 0
    convolutionPairPeaks_noRedundancy = {}
    pair_groups = {}
    for pair in convolutionPairPeaks_list:
        new_group = True
        if current_group == 0:
            current_group += 1
            minmax_x = [pair[0], pair[0]]
            minmax_y = [pair[1], pair[1]]
            convolutionPairPeaks_noRedundancy[current_group] = [pair]
            pair_groups[pair] = current_group
            continue
        adduct = 28
        if minmax_x[0] - adduct < pair[0] < minmax_x[1] + adduct:
            if minmax_y[0] - adduct < pair[1] < minmax_y[1] + adduct:
                new_group = False
        if not new_group:
            if pair[0] < minmax_x[0]:
                minmax_x[0] = pair[0]
            if pair[0] > minmax_x[1]:
                minmax_x[1] = pair[0]

            if pair[1] < minmax_y[0]:
                minmax_y[0] = pair[1]
            if pair[1] > minmax_y[1]:
                minmax_y[1] = pair[1]
            convolutionPairPeaks_noRedundancy[current_group].append(pair)
            pair_groups[pair] = current_group
        else:
            current_group += 1
            minmax_x = [pair[0], pair[0]]
            minmax_y = [pair[1], pair[1]]
            convolutionPairPeaks_noRedundancy[current_group] = [pair]
            pair_groups[pair] = current_group
            continue
    min_element = offset
    max_element = offset
    alldistances = {}
    final_group = {}
    for pair in convolutionPairPeaks_list:
        if round(abs(pair[1] - pair[0]), 3) not in alldistances:
            alldistances[round(abs(pair[1] - pair[0]), 3)] = [pair]
        else:
            alldistances[round(abs(pair[1] - pair[0]), 3)].append(pair)
    # if  abs(round(abs(pair[1]-pair[0]),3) - offset) < e:
    # 	final_group.append(round(abs(pair[1]-pair[0]),3))
    sorted_alldistances = sorted(alldistances)
    last_distance = 0
    aa_group = False
    last_group_status = False
    final_group_distanes = []
    for i in range(len(sorted_alldistances)):
        if i == len(sorted_alldistances) - 1:
            if last_group_status:
                if abs(distance - offset) > e:
                    break
                else:
                    final_group_distanes.append(distance)
                    break
        distance = sorted_alldistances[i]
        if abs(distance - offset) < e:
            aa_group = True
        if abs(distance - last_distance) > e:
            if last_group_status:
                if abs(distance - offset) > e:
                    break
            final_group_distanes = [distance]

        else:
            final_group_distanes.append(distance)
        last_distance = distance
        last_group_status = aa_group
    final_groups = {}

    for distance in final_group_distanes:
        final_groups[distance]
        for pair in alldistances[distance]:
            final_groups[pair_groups[pair]] = convolutionPairPeaks_noRedundancy[pair_groups[pair]]

    return final_groups
    '''
