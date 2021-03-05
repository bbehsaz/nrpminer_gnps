import os
import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import re

def read_gbk_file_get_features(gb_file):
    from Bio import SeqIO
#     gb_file = "LOIC01000001.1.cluster017.gbk"
    for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank") :
        # now do something with the record
#         print "Name %s, %i features" % (gb_record.name, len(gb_record.features))
#         print repr(gb_record.seq)
        return gb_record



def index_genbank_features(gb_record, feature_type, qualifier) :
    import os
    import Bio
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio import SeqIO
    import re
    answer = dict()
    for (index, feature) in enumerate(gb_record.features) :
        if feature.type==feature_type :
            if qualifier in feature.qualifiers :
                #There should only be one locus_tag per feature, but there
                #are usually several db_xref entries
                for value in feature.qualifiers[qualifier] :
                    if value in answer :
                        print "WARNING - Duplicate key %s for %s features %i and %i" \
                           % (value, feature_type, answer[value], index)
                    else :
                        answer[value] = index
    return answer
def index_genbank_clusters_clusterID(gb_record, feature_type, qualifier) :
    import os
    import Bio
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio import SeqIO
    import re
    answer = dict()
    for (index, feature) in enumerate(gb_record.features) :
        if feature.type==feature_type :
            print feature
            if qualifier in feature.qualifiers :
                print feature.qualifiers['note'] 
                clusterID = feature.qualifiers['note'][0]
                print clusterID
                #There should only be one locus_tag per feature, but there
                #are usually several db_xref entries
                answer[clusterID] = index
#                 for value in feature.qualifiers[qualifier] :
#                     if value in answer :
#                         print "WARNING - Duplicate key %s for %s features %i and %i" \
#                            % (value, feature_type, answer[value], index)
#                     else :
#                         answer[value] = index
    return answer
def index_nrps_genbank_features(gb_record, feature_type, qualifier) :
    import os
    import Bio
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio import SeqIO
    import re
    answer = dict()
    predictions = dict()
    for (index, feature) in enumerate(gb_record.features) :
        if feature.type==feature_type :
            if qualifier in feature.qualifiers :
                #There should only be one locus_tag per feature, but there
                #are usually several db_xref entries
                if "Type" in feature.qualifiers['sec_met']:
                        product_type = info.split(":")[1][1:]
                        if product_type !='nrps':
                            continue
                for value in feature.qualifiers['locus_tag'] :
                    if value in answer :
                        print "WARNING - Duplicate key %s for %s features %i and %i" \
                           % (value, feature_type, answer[value], index)
                    else :
                        answer[value] = index
                        predictions[value] = feature.qualifiers['aSProdPred'][0]
    return answer, predictions
