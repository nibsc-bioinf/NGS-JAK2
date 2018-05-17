import pysam
import os
import csv
import subprocess
import re

'''
NB: cigartuple codes: 0 = match, 1 = ins, 2 = del, 3 = ref_skip, 4 = soft clip,
    5 = hard clip, 6 = pad, 7 = equal, 8 = diff, 9 = back.
'''


class MutationReadContext:
    '''
    Class for creating objects to store information about the read context.
    i.e. the position of the mutation site from various read features.
    '''
	
    def __init__(self, project_id, percent, run, read_id, orientation, allele, 
                start, hc_start, sc_start, sc_end, hc_end, end, ins_count, 
                closest_ins_dist, del_count, closest_del_dist, mut_count, 
                snp_count):
        self.project_id = project_id
        self.percent = percent
        self.run = run
        self.read_id = read_id
        self.orientation = orientation
        self.allele = allele
        self.start = start
        self.hc_start = hc_start
        self.sc_start = sc_start
        self.sc_end = sc_end
        self.hc_end = hc_end
        self.end = end
        self.ins_count = ins_count
        self.closest_ins_dist = closest_ins_dist
        self.del_count = del_count
        self.closest_del_dist = closest_del_dist
        self.mut_count = mut_count
        self.snp_count = snp_count
				 
			
class MutationContextParser:
    '''
    Parses all BAM files in the input folder and generates a TSV file for the
    mutation context of all reads which map to the region of interest 
    (position 1226) 
    '''
 	
    def get_reads_of_interest_from_BAM(self):    
        # finds all BAM files in input folder and outputs all the reads of 
        # interest to a new BAM file in the output folder 					
        for filename in os.listdir("data/input/"):
            if filename.endswith(".bam"):
                # indexes BAM files waiting to be processed
                pysam.sort("data/input/" + filename)
                pysam.index("data/input/" + filename)
                samfile = pysam.AlignmentFile("data/input/" + filename, "rb")
                reads_of_interest = pysam.AlignmentFile("data/output/" + 
                    filename + "_ROI.bam", "wb", template=samfile)
                # fetches all read which include postion 1226
                for read in samfile.fetch("JAK2", 1225, 1226): 
                    reads_of_interest.write(read)
                reads_of_interest.close()
                samfile.close()
                
    def remove_duplicates(self):
        # removes duplicated reads using Picard tools MarkDuplicates
        for filename in os.listdir("data/output/"):
                if filename.endswith(".bam"):
                    cmd = "java -jar /usr/local/share/picard/picard.jar " + \
                            "MarkDuplicates REMOVE_DUPLICATES=true " + \
                            "I=data/output/" + filename + " O=data/output/" + \
                            filename[:-3] + "rm.bam METRICS_FILE=samp.dups " + \
                            "ASSUME_SORTED=true"
                    print(cmd)
                    run = subprocess.Popen(cmd, shell=True)
                    run.wait()
                            
    def sort_and_index_BAMs(self):
        # sorts and indexes all BAM files in the output folder            
        for filename in os.listdir("data/output/"):
            if filename.endswith(".rm.bam"):
                pysam.sort("data/output/" + filename)
                pysam.index("data/output/" + filename)
                
    def get_orientation(self, read):
        # determines whether the read is forward or reverse
        if read.is_read1:
            orientation = 'F'
        else:
            orientation = 'R'
        return orientation
                
    def get_allele(self, read, cigartuples):
        # checks if read starts with soft clipping
        if cigartuples[0][0] == 4: 
            # finds allele at position 1226 offset by amount of soft clipping
            try:
                allele = read.seq[read.get_reference_positions().index(1225) +  
                            cigartuples[0][1]]
            except:
                allele = "del"
        # checks if read starts with hard clipping followed by soft clipping
        elif cigartuples[0][0] == 5 and cigartuples[1][0] == 4:
            # finds allele at position 1226 offset by soft clipping
            try:
                allele = read.seq[read.get_reference_positions().index(1225) +  
                            cigartuples[1][1]]
            except:
                allele = "del"
        else:
            # else get the allele at position 1226
            try:
                allele = read.seq[read.get_reference_positions().index(1225)]
            except:
                allele = "del"
        return allele
        
    def update_allele(self, read, cigartuples, mut_allele_offset, allele):
        pos = 1225 + mut_allele_offset
        # checks if read starts with soft clipping
        if cigartuples[0][0] == 4: 
            # finds allele at position 1226 offset by amount of soft clipping
            try:
                allele = read.seq[read.get_reference_positions().index(pos) +  
                            cigartuples[0][1]]
            except:
                allele = 'G'
                print(read.query_name, mut_allele_offset)
        # checks if read starts with hard clipping followed by soft clipping
        elif cigartuples[0][0] == 5 and cigartuples[1][0] == 4:
            # finds allele at position 1226 offset by soft clipping
            try:
                allele = read.seq[read.get_reference_positions().index(pos) +  
                            cigartuples[1][1]]
            except:
                allele = 'G'
                print(read.query_name, mut_allele_offset)
        else:
            # else get the allele at position 1226
            try:
                allele = read.seq[read.get_reference_positions().index(pos)]
            except:
                allele = 'G'
                print(read.query_name, mut_allele_offset)
        return allele
    
    def get_position_from_start(self, read, cigartuples):
        # checks whether read starts with soft clipping        
        if cigartuples[0][0] == 4: 
            # sets postion of first nucleotide in sequence offset by target 
            # nucleotide and soft clipping
            position = 1225 - read.get_reference_positions()[0] + \
                        cigartuples[0][1]
        # checks if read starts with hard clipping followed by soft clipping
        elif cigartuples[0][0] == 5 and cigartuples[1][0] == 4:
            # finds allele at position 1226 offset by soft clipping
            position = 1225 - read.get_reference_positions()[0] + \
                        cigartuples[1][1]
        else:
            # else sets postion of first nucleotide in sequence offset by 
            # target nucleotide
            position = 1225 - read.get_reference_positions()[0]
        return position
        
    def get_position_from_sc_start(self, start, cigartuples):
        sc_start = start
        # checks whether read starts with soft clipping
        if cigartuples[0][0] == 4:
            # offset sc_start by amount of soft clipping
            sc_start -= cigartuples[0][1]
        # checks if read starts with hard clipping followed by soft clipping
        elif cigartuples[0][0] == 5 and cigartuples[1][0] == 4:
            # offset sc_start by amount of soft clipping
            sc_start -= cigartuples[1][1]
        return sc_start
    
    def get_position_from_hc_start(self, start, cigartuples):
        hc_start = start
        # checks whether read starts with hard clipping
        if cigartuples[0][0] == 5:
            # start position changed to start of read before hard clipping
            start += cigartuples[0][1]
        return hc_start, start
        
    def get_position_from_end(self, read, cigartuples):
        # checks whether read ends with soft clipping    
        if cigartuples[-1][0] == 4: 
            # sets postion of last nucleotide in sequence offset by target 
            # nucleotide and soft clipping
            position = read.get_reference_positions()[-1] - 1225 + \
                        cigartuples[-1][1]
        # checks if read ends with hard clipping followed by soft clipping
        elif cigartuples[-1][0] == 5 and cigartuples[-2][0] == 4:
            # finds allele at position 1226 offset by soft clipping
            position = read.get_reference_positions()[-1] - 1225 + \
                        cigartuples[-2][1]
        else:
            # else sets postion of last nucleotide in sequence offset by 
            # target nucleotide
            position = read.get_reference_positions()[-1] - 1225
        return position
    
    def get_position_from_sc_end(self, end, cigartuples):
        sc_end = end
        # checks whether read ends with soft clipping
        if cigartuples[-1][0] == 4:
            # offset sc_end by amount of soft clipping
            sc_end -= cigartuples[-1][1]
        # checks if read ends with hard clipping followed by soft clipping
        elif cigartuples[-1][0] == 5 and cigartuples[-2][0] == 4:
            # offset sc_end by amount of soft clipping
            sc_end -= cigartuples[-2][1]
        return sc_end
        
    def get_position_from_hc_end(self, end, cigartuples):
        hc_end = end
        # checks whether read ends with hard clipping
        if cigartuples[-1][0] == 5:
            # end position changed to end of read before hard clipping
            end += cigartuples[-1][1]
        return hc_end, end
        
    def get_ins_count(self, cigartuples):
        ins_count = 0
        # loops through all cigartuples
        for cigartuple in cigartuples:
            # if cigartuple is an insert, adds 1 to the insert count
            if cigartuple[0] == 1:
                ins_count += 1
        return ins_count
        
    def get_ins_dist(self, read, cigartuples):
        # list for storing insert distances from the read start
        insert_dists_from_start = []
        # list for storing insert distances from the mutation site
        insert_dists_from_target = []
        # insert distance from start initialised to zero
        ins_dist_from_start = 0
        # closest insert distance initialised to a large number
        closest_ins_dist = 99999
        for cigartuple in cigartuples:
            # if the cigar tuple is not an insert
            if cigartuple[0] != 1:
                # increments distance from start by length of feature
                ins_dist_from_start += cigartuple[1]
            # if the cigar tuple is an insert
            if cigartuple[0] == 1:
                # set the insert distance from start
                ins_dist = ins_dist_from_start
                # add the insert to list with size of insert (cigartuple[1])
                insert_dists_from_start.append((ins_dist, cigartuple[1]))
                # increments distance from start by length of insert
                ins_dist_from_start += cigartuple[1]
        # for each insert
        for insert in insert_dists_from_start:
            # calculate the distance from the mutation site
            # if read starts with soft clipping
            if cigartuples[0][0] == 4:
                # offset insert position by soft clipping length
                dist_from_target = read.get_reference_positions()[
                                    insert[0] - cigartuples[0][1]] - 1225
            elif cigartuples[0][0] == 5 and cigartuples[1][0] == 4:
                # offset insert position by soft clipping length
                dist_from_target = read.get_reference_positions()[
                                    insert[0] - cigartuples[1][1]] - 1225
            elif cigartuples[0][0] == 5:
                # offset insertion position by hard clipping length
                dist_from_target = read.get_reference_positions()[
                                    insert[0] - cigartuples[0][1]] - 1225
            else:
                dist_from_target = read.get_reference_positions()[
                                    insert[0]] - 1225
            # count for total size of all inserts prior to mutation site
            mut_allele_offset = 0
            # if the insert appears before the mutation site
            if dist_from_target < 0:
                # change the distance so that it is the end of the insert
                dist_from_target += insert[1] - 1 
                mut_allele_offset += insert[1]
            # add the absolute distance from mutation site to the list
            insert_dists_from_target.append(abs(dist_from_target))
        # for each item in the list
        for ins_dist in insert_dists_from_target:
            # compare the distance to the current closest distance
            if ins_dist < closest_ins_dist:
                # if it is closer than current closest, update current closest
                closest_ins_dist = ins_dist
        return closest_ins_dist, mut_allele_offset
        
    def get_del_count(self, cigartuples):
        del_count = 0
        # loops through all cigartuples
        for cigartuple in cigartuples:
            # if cigartuple is a deletion, adds 1 to the deletion count
            if cigartuple[0] == 2:
                del_count +=1
        return del_count
        
    def get_del_dist(self, read, cigartuples):
        # list for storing deletion distances from the read start
        deletion_dists_from_start = []
        # list for storing deletion distances from the mutation site
        deletion_dists_from_target = []
        # deletion distance from start initialised to zero
        del_dist_from_start = 0
        # closest deletion distance initialised to a large number
        closest_del_dist = 99999
        for cigartuple in cigartuples:
            # if the cigar tuple is not a deletion
            if cigartuple[0] != 2:
                # increments distance from start by length of feature
                del_dist_from_start += cigartuple[1]
            # if the cigar tuple is a deletion
            if cigartuple[0] == 2:
                # set the deletion distance from start
                del_dist = del_dist_from_start
                # add the deletion to list with size (cigartuple[1])
                deletion_dists_from_start.append((del_dist, cigartuple[1]))
                # increments distance from start by length of deletion
                del_dist_from_start += cigartuple[1]
        # for each deletion
        for deletion in deletion_dists_from_start:
            # calculate the distance from the mutation site
            # if read starts with soft clipping
            if cigartuples[0][0] == 4:
                # offset deletion position by soft clipping length
                dist_from_target = read.get_reference_positions()[
                                    deletion[0] - cigartuples[0][1]] - 1225
            elif cigartuples[0][0] == 5 and cigartuples[1][0] == 4:
                # offset deletion position by soft and hard clipping length
                dist_from_target = read.get_reference_positions()[
                                    deletion[0] - cigartuples[0][1] - 
                                    cigartuples[1][1]] - 1225
            elif cigartuples[0][0] == 5:
                # offset deletion position by hard clipping length
                dist_from_target = read.get_reference_positions()[
                                    deletion[0] - cigartuples[0][1]] - 1225
            else:
                dist_from_target = read.get_reference_positions()[
                                        deletion[0]] - 1225
            # if the deletion appears before the mutation site
            if dist_from_target < 0:
                # change the distance so that it is the end of the deletion
                dist_from_target += deletion[1] - 1 
            # add the absolute distance from mutation site to the list
            deletion_dists_from_target.append(abs(dist_from_target))
        # for each item in the list
        for del_dist in deletion_dists_from_target:
            # compare the distance to the current closest distance
            if del_dist < closest_del_dist:
                # if it is closer than current closest, update current closest
                closest_del_dist = del_dist
        return closest_del_dist
    
    def get_mut_count(self, read, allele, ins_count, del_count):
        count = 0
        # gets MD tag and splits out into numbers and letters
        md_tag = (re.split('(\d+)', read.get_tag("MD")))
        for section in md_tag:
            # if part of MD tag is a SNP (ignores deletions which start with ^)
            if section.isalpha():
                count += 1
        # # if read contains mutant allele, minus 1 from SNP count total
        count += ins_count + del_count
        if allele == "T":
            count -= 1
        return count
        
    def get_snp_count(self, read, allele):
        count = 0
        # gets MD tag and splits out into numbers and letters
        md_tag = (re.split('(\d+)', read.get_tag("MD")))
        for section in md_tag:
            # if part of MD tag is a SNP (ignores deletions which start with ^)
            if section.isalpha():
                count += 1
        # if read contains mutant allele, minus 1 from SNP count total
        if allele == "T":
            count -= 1
        return count    
        
    def write_to_tsv(self, mutation_context_list):
        with open ("data/output/tsv_files/results.tsv", 'w') as tsvfile:
            f = csv.writer(tsvfile, delimiter='\t')
            # writes each of the mutation context features to the file
            for mc in mutation_context_list:
                f.writerow([mc.project_id, mc.percent, mc.run, mc.read_id, 
                            mc.orientation, mc.start, mc.hc_start, mc.sc_start,
                            mc.sc_end, mc.hc_end, mc.end, mc.allele, 
                            mc.ins_count, mc.closest_ins_dist, mc.del_count,
                            mc.closest_del_dist, mc.mut_count, mc.snp_count])
        print("Results written to disk")
        
    def parse_BAM_files(self):
        self.get_reads_of_interest_from_BAM()
        self.remove_duplicates()
        self.sort_and_index_BAMs()
        # list for storing the mutation context for each read
        mutation_context_list = []
        # finds all BAM files in output folder        
        for filename in os.listdir("data/output/"):
            if filename.endswith(".rm.bam"):
                project_id = filename[:3]
                percent = filename.split('-', 1)[1].split('per', 1)[0]
                run = filename.split('rep', 1)[1].split('.', 1)[0]
                samfile = pysam.AlignmentFile("data/output/" + filename, "rb")
                for read in samfile:
                    # gets the cigar tuples which describe the cigar string
                    cigartuples = read.cigartuples
                    orientation = self.get_orientation(read)
                    allele = self.get_allele(read, cigartuples)
                    start = self.get_position_from_start(read, cigartuples)
                    sc_start = self.get_position_from_sc_start(start,
                                                                cigartuples)
                    hc_start, start = self.get_position_from_hc_start(start, 
                                                                cigartuples)
                    end = self.get_position_from_end(read, cigartuples)
                    sc_end = self.get_position_from_sc_end(end, cigartuples)
                    hc_end, end = self.get_position_from_hc_end(end, 
                                                                cigartuples)
                    ins_count = self.get_ins_count(cigartuples)
                    # sets closest_ins_dist to - for reads with no inserts
                    closest_ins_dist = '-'
                    if ins_count > 0:
                        closest_ins_dist, mut_allele_offset = \
                                           self.get_ins_dist(read, cigartuples)
                        
                        allele = self.update_allele(read, cigartuples, 
                                                    mut_allele_offset, allele)
                    del_count = self.get_del_count(cigartuples)
                    # sets closest_del_dist to - for reads with no dels
                    closest_del_dist = '-'
                    if del_count > 0:
                        closest_del_dist = self.get_del_dist(read, cigartuples)
                    mut_count = self.get_mut_count(read, allele, ins_count,
                                                    del_count)
                    snp_count = self.get_snp_count(read, allele)
                    # creates a MutationReadContext object for each read
                    mutation_context = MutationReadContext(project_id, percent, 
                                        run, read.query_name, orientation, 
                                        allele, start, hc_start, sc_start, 
                                        sc_end, hc_end, end, ins_count, 
                                        closest_ins_dist, del_count, 
                                        closest_del_dist, mut_count, snp_count)
                    mutation_context_list.append(mutation_context)
                #self.write_to_tsv(mutation_context_list, filename)
                samfile.close()
        self.write_to_tsv(mutation_context_list)
        print("All files processed")

# runs the script        
if __name__ == "__main__":
    parser = MutationContextParser()
    parser.parse_BAM_files()