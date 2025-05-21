#include <htslib/sam.h>
#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cstring>
#include <cstdlib>
#include <tuple>
#include <math.h>

using namespace std;
struct BedInterval {
    std::string chrom;
    int start;
    int end;
    std::string label() const {
        return chrom + ":" + std::to_string(start) + "-" + std::to_string(end);
    }
};

std::unordered_map<std::string, std::vector<BedInterval>> read_bed_file(const std::string& bed_path, vector<string> &labels) {
    std::unordered_map<std::string, std::vector<BedInterval>> bed_map;
    std::ifstream bed(bed_path);
    std::string line;

    while (std::getline(bed, line)) {
        std::istringstream iss(line);
        BedInterval bi;
        iss >> bi.chrom >> bi.start >> bi.end;
	labels.push_back(bi.label());
        bed_map[bi.chrom].push_back(bi);
    }

    for (auto& [chrom, vec] : bed_map) {
        std::sort(vec.begin(), vec.end(), [](const BedInterval& a, const BedInterval& b) {
            return a.start < b.start;
        });
    }

    return bed_map;
}

std::vector<BedInterval> find_overlapping_intervals(const std::vector<BedInterval>& intervals, int read_start, int read_end) {
    std::vector<BedInterval> results;
    auto it = std::lower_bound(intervals.begin(), intervals.end(), read_start,
        [](const BedInterval& a, int val) {
            return a.end <= val;
        });

    while (it != intervals.end() && it->start < read_end) {
        results.push_back(*it);
        ++it;
    }

    return results;
}

string GetMMDeltaStr(const std::string & mm_tag) {
    size_t first_comma = mm_tag.find(',');
    if (first_comma == std::string::npos) return "";

    std::string delta_str = mm_tag.substr(first_comma + 1);
    return delta_str;
}

void GetMMDeltas(string delta_str, vector<int> &deltas) {
  if (delta_str.empty()) { deltas.clear(); return;}
    stringstream delta_strm(delta_str);
    int pos=0;
    int strPos=0;
    while(delta_strm) {
      int delta=0;
      char c;
      if ( !(delta_strm >> delta) ) {
	break;
      }
      deltas.push_back(delta);
      delta_strm.get();
    }  
    return;
}

void get_mm_annot(const std::string& mm_tag, const std::string& seq, char base, string &annot) {

  string deltaStr = GetMMDeltaStr(mm_tag);
  vector<int> deltas;
  GetMMDeltas(deltaStr, deltas);
  int deltaIndex=0;
  int skip=0;
  char annotChar;
  for (auto i =0 ; i < seq.size() and deltaIndex < deltas.size(); i++) {
    
    if (seq[i] == base) {
      if (skip == deltas[deltaIndex]) {
	annot.push_back('*');
	skip=0;
	deltaIndex++;
      }
      else {
	skip++;
      }
    }
    annot.push_back(seq[i]);
  }
}


std::vector<int> parse_mm_tag(const std::string& mm_tag, const std::string& seq, char base, int strand) {
    std::vector<int> mod_positions;
    //    if (mm_tag.empty() || mm_tag[0] != base) return mod_positions;
    string mmDeltaStr = GetMMDeltaStr(mm_tag);
    if (mmDeltaStr.empty())  return mod_positions;
    int delta = 0;
    size_t i = 0;
    stringstream delta_strm(mmDeltaStr);
    int pos=0;
    int strPos=0;
    while(delta_strm) {
      int delta=0;
      char c;
      if ( !(delta_strm >> delta) ) {
	break;
      }
      int nc=0;
      while (strPos < seq.size() and nc < delta) {
	if (seq[strPos] == base) { nc++;}
	strPos++;
      }
      while (strPos < seq.size() and seq[strPos] != base) { strPos++;}

      
      pos+=delta;
      if (delta != 0) {
	if (strand == 0) {
	  mod_positions.push_back(strPos);
	}
	else {
	  mod_positions.push_back(seq.size() - strPos - 1);
	}
      }
      strPos++;
      delta_strm.get();
    }

    return mod_positions;
}

class Interval {
public:
  string chrom;
  int start;
  int end;
};

bool NegComp(const int &a, const int &b) {
  return b < a;
}

int get_hp_tag(bam1_t* b, bool& found) {
    found = false;
    uint8_t* hp_data = bam_aux_get(b, "HP");
    if (!hp_data) return -1;

    char type = hp_data[0];
    if (type == 'i') {
        found = true;
        return bam_aux2i(hp_data);
    } else if (type == 'I') {
        found = true;
        return static_cast<int>(bam_aux2i(hp_data));
    } else {
        // Unexpected type
        return -1;
    }
}

void GetHapIndices(int hap, vector<int> &hapIdx) {
  hapIdx.clear();
  if (hap == 1) {
    hapIdx.push_back(0);
  }
  else if (hap == 2) {
    hapIdx.push_back(1);
  }
  else {
    hapIdx.push_back(0);
    hapIdx.push_back(1);
  }
}

void SummaryStats(vector<pair<int,int> > &vals, int startOffset, int endOffset, float &mean, float &sd) {
  float sum=0,sumsq=0;
  if (vals.size() <= endOffset) {
    mean=sd=0;
    return;
  }
  for (int i=startOffset; i < vals.size() - endOffset; i++) {
    int v=vals[i].first;
    sum+=v; sumsq += v*v;
  }
  int n=(vals.size() - startOffset - endOffset);
  mean=((float)sum)/n;
  sd = sqrt(sumsq/n - mean*mean);  
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: countme <input.bam> <input.bed>\n";
        return 1;
    }
    vector<string> labels;
    auto bed_intervals = read_bed_file(argv[2], labels);

    // map: interval label -> (sum of per-read avg methylation, count of reads contributing)
    vector<std::unordered_map<std::string, std::pair<int, int>> > interval_sums(2);
    vector<std::unordered_map<std::string, std::pair<int, int>> > interval_lengths(2);
    vector<std::unordered_map<std::string, vector<int> > > interval_length_vect(2);
    vector<std::unordered_map<std::string, vector<int> > > interval_cpgCount(2);
    vector<std::unordered_map<std::string, vector<int> > > interval_cpgMeth(2);    
    std::unordered_map<std::string, BedInterval> interval_bed;
    
    samFile* in = sam_open(argv[1], "r");
    if (!in) {
        std::cerr << "Error opening BAM file\n";
        return 1;
    }

    bam_hdr_t* hdr = sam_hdr_read(in);
    bam1_t* b = bam_init1();
    int proc_read=0;

    cerr << "Running countme v0.61" << endl;
    // Initialie counters for all intervals
    vector<int> hi;
    hi.push_back(0);
    hi.push_back(1);
    for (auto& chrom: bed_intervals) {
      for (auto& bed: chrom.second) {
	string label = bed.label();
	interval_bed[label] = bed;
	for (auto h: hi) {
	  interval_sums[h][bed.label()].first =0;
	  interval_sums[h][bed.label()].second = 0;
	  interval_lengths[h][bed.label()].first = 0;
	  interval_lengths[h][bed.label()].second = 0;
	  interval_length_vect[h][bed.label()] = vector<int>();
	  interval_cpgCount[h][bed.label()] = vector<int>();
	  interval_cpgMeth[h][bed.label()] = vector<int>();	  
	}
      }
    }
    while (sam_read1(in, hdr, b) >= 0) {
      
        if (b->core.flag & BAM_FUNMAP) continue;
	++proc_read;
	if (proc_read % 100000 == 0) {
	  cerr << "Proc " << proc_read << endl;
	}
        std::string chrom = hdr->target_name[b->core.tid];
        int ref_start = b->core.pos; 
        int ref_end = bam_endpos(b);
	
        std::string read_name = bam_get_qname(b);
	bool hpFound;
	int haplotype = get_hp_tag(b, hpFound);
	vector<int> hapIndices;
	GetHapIndices(haplotype, hapIndices);
        auto chrom_it = bed_intervals.find(chrom);
        if (chrom_it == bed_intervals.end()) continue;

        uint8_t* mm_data = bam_aux_get(b, "MM");
        uint8_t* ml_data = bam_aux_get(b, "ML");
	if (mm_data == NULL) {
	  mm_data = bam_aux_get(b, "Mm");
	}
	if (ml_data == NULL) {
	  ml_data = bam_aux_get(b, "Ml");
	}
	
        if (!mm_data || !ml_data) continue;

        std::string mm_str = bam_aux2Z(mm_data);
        int read_len = b->core.l_qseq;
	if (read_len == 0) { continue;}
        std::vector<uint8_t> meth_status(read_len, 0);

        std::string seq(read_len, 'N');
        uint8_t* seq_ptr = bam_get_seq(b);
        for (int i = 0; i < read_len; ++i)
            seq[i] = seq_nt16_str[bam_seqi(seq_ptr, i)];

        uint8_t type = ml_data[0];
        if (type != 'B' || ml_data[1] != 'C') continue;
        uint32_t count;
        std::memcpy(&count, ml_data + 2, sizeof(uint32_t));
	vector<int> ml_array( count);
	for (auto i=0; i < count; i++ ){ ml_array[i] = ml_data[i+6];}
	  //        const uint8_t* ml_array = reinterpret_cast<uint8_t*>(ml_data + 6);

	char base = 'C';
	int flag = b->core.flag;
	int strand =0;
	if (flag & BAM_FREVERSE) {
	  strand = 1;
	}

	if (strand == 1 ) {
	  string seq_rc(seq);
	  int n=seq.size();
	  for (auto i = 0; i < seq.size(); i++) {
	    if (seq[i] == 'A') { seq_rc[n-i-1] = 'T';}
	    if (seq[i] == 'T') { seq_rc[n-i-1] = 'A';}
	    if (seq[i] == 'C') { seq_rc[n-i-1] = 'G';}
	    if (seq[i] == 'G') { seq_rc[n-i-1] = 'C';}	    
	  }
	  seq=seq_rc;
	}

        auto pos_c = parse_mm_tag(mm_str, seq, base, strand);

        std::vector<int> read_to_ref(read_len, -1);
        uint32_t* cigar = bam_get_cigar(b);
        int read_pos = 0, ref_pos = b->core.pos;

        for (uint32_t i = 0; i < b->core.n_cigar; ++i) {
            uint32_t op = bam_cigar_op(cigar[i]);
            uint32_t len = bam_cigar_oplen(cigar[i]);

            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                for (uint32_t j = 0; j < len; ++j) {
		  if (strand == false) {
                    if (read_pos < read_len) read_to_ref[read_pos++] = ref_pos++;
		  }
		  else {
		    ++read_pos;
                    if (read_pos > 0) read_to_ref[seq.size() - read_pos] = ref_pos++;
		  }
                }
            } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
                ref_pos += len;
            } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
	      for (uint32_t j=0; j < len; ++j) {
		  if (strand == false) {
                    if (read_pos < read_len) read_to_ref[read_pos++] = ref_pos;
		  }
		  else {
		    ++read_pos;
                    if (read_pos > 0) read_to_ref[seq.size() - read_pos] = ref_pos;
		  }
	      }		
            }
        }

	vector<int>::iterator readAlnStart = read_to_ref.begin(),  readAlnEnd = read_to_ref.end();
	while(readAlnStart != read_to_ref.end() and *readAlnStart == -1) { readAlnStart++;}
	while(readAlnEnd != readAlnStart ) {
	  vector<int>::iterator tmpIt = readAlnEnd;
	  tmpIt--;
	  if (*tmpIt != -1) break;
	  readAlnEnd--;
	}
        size_t m_idx = 0;
	
        for (int pos : pos_c) {
            if (m_idx >= count) break;
	    assert(pos < meth_status.size());
            meth_status[pos] = ml_array[m_idx++];
        }


        auto overlapping = find_overlapping_intervals(chrom_it->second, ref_start, ref_end);
	int nMatched=0;
        for (const auto& bed : overlapping) {
            std::vector<int> relevant;
	    int methylated = 0;
	    int total_c = 0;
	    int minOverhang=100;
	    vector<int>::iterator readIntvStart, readIntvEnd;
	    if (ref_start + minOverhang >  bed.start) {
	      //	      cout << read_name << " starts within bed interval" <<  bed.end - bed.start << endl;
	      continue;
	    }
	    if (ref_end -minOverhang< bed.end) {
	      //	      cout << read_name << " ends before end of bed interval " <<  bed.end - bed.start << endl;
	      continue;
	    }
	    
	    if (strand == 0) {
	      readIntvStart = lower_bound(readAlnStart, readAlnEnd, bed.start-10);
	      readIntvEnd   = upper_bound(readAlnStart, readAlnEnd, bed.end+10);
	    }
	    else {
	      readIntvStart = upper_bound(readAlnStart, readAlnEnd, bed.end+10, NegComp);
	      readIntvEnd   = lower_bound(readAlnStart, readAlnEnd, bed.start-10, NegComp);
	    }
	    if (readIntvStart == read_to_ref.end() or readIntvEnd == read_to_ref.end() ) { continue;}

	    if ((strand == 0 and *readIntvStart > bed.start+10) or (strand == 1 and *readIntvEnd > bed.start+10)) {
	      //	      cout << "Read starting inside interval" << endl;
	      continue;
	    }
	    int readIntvStartIdx = readIntvStart - read_to_ref.begin()+10;
	    int readIntvEndIdx = readIntvEnd - read_to_ref.begin()-10;
	    if (readIntvStartIdx >= readIntvEndIdx) {
	    }

	    else {
	      nMatched++;
	    }
            for (int i = readIntvStartIdx; i < readIntvEndIdx; ++i) {
                int rp = read_to_ref[i];
                if (rp >= bed.start && rp < bed.end && i < seq.size()-1 and
		    (( strand == 0 and seq[i] == 'C' and seq[i+1] == 'G')  or
		     ( strand == 1 and i > 0 and seq[i-1] == 'C' and seq[i] == 'G'))) {
		  total_c++;
                    if (meth_status[i] >= 127) {
		      ++methylated;
                    }
                }
            }

	    for (auto h: hapIndices) {
	      interval_sums[h][bed.label()].first += methylated;
	      interval_sums[h][bed.label()].second += total_c;
	      //	      cout << read_name<< " hap: " << h << " adding length " << readIntvEndIdx - readIntvStartIdx<< " " << readIntvEndIdx  << " " << readIntvStartIdx << endl;
	      interval_lengths[h][bed.label()].first += readIntvEndIdx - readIntvStartIdx;
	      interval_cpgCount[h][bed.label()].push_back(total_c);
	      interval_cpgMeth[h][bed.label()].push_back(methylated);
	      interval_length_vect[h][bed.label()].push_back(readIntvEndIdx - readIntvStartIdx);
	      interval_lengths[h][bed.label()].second += 1;
	    }
	    interval_bed[bed.label()] = bed;
	}
	//	cout << "For " << read_name << " ovp " << overlapping.size() << " matched: " << nMatched << endl;
    }

    for (auto label: labels) {
      vector<double> avgMethRatio(2,0), avgCpGCount(2,0), avgLen(2,0), nReads(2,0);
      bool didDelete=false;
      for (auto i =0; i < 2; i++) {
	if (interval_sums[i].find(label) != interval_sums[i].end()) {
	  const auto& [meSum, cpgCount] = interval_sums[i][label];
	  const auto& [lenSum, readCount] = interval_lengths[i][label];
	  double meanMeth = cpgCount == 0 ? -1.0 : ((float)meSum) / cpgCount;
	  double meanLen = readCount == 0 ? -1.0 : ((float) lenSum) / readCount;
	  nReads[i] = readCount;
	  avgMethRatio[i] = meanMeth;
	  avgCpGCount[i] = cpgCount;
	  avgLen[i] = meanLen;
	}
	else {
	  nReads[i] = 0;
	  avgMethRatio[i] = 0;
	  avgCpGCount[i] = 0;
	  avgLen[i] = 0;
	}
	//
	// Check lengths for suspect values.
	//
	vector<int> toDelete;	
	if (interval_length_vect[i][label].size() > 10) {
	  vector<pair<int,int> >lengths;
	  for (int li =0; li < interval_length_vect[i][label].size(); li++) {
	    lengths.push_back(pair<int,int>(interval_length_vect[i][label][li], li));
	  }
	  std::sort(lengths.begin(), lengths.end());

	  float loMean, loSD, fullMean, fullSD, hoMean, hoSD;
	  SummaryStats(lengths, 1, 0, loMean, loSD);
	  SummaryStats(lengths, 0, 1, hoMean, hoSD);
	  int nValue=lengths.size();
	  for (int l=0; l < nValue; l++) {
	    if (interval_length_vect[i][label][l] > hoMean + 3*hoSD + 100 or
		interval_length_vect[i][label][l] < hoMean - 3*hoSD - 100) {		
	      //	      cerr << "DELETING: " << interval_length_vect[i][label][l] << endl;
	      toDelete.push_back(l);
	    }
	  }

	  int cur=0;
	  for (int j=0; j < interval_length_vect[i][label].size(); j++) {
	    bool deleteThis=false;
	    for (auto d: toDelete) {
	      if (j == d) {
		deleteThis=true;
		break;
	      }
	    }
	    if (deleteThis == false) {
	      interval_length_vect[i][label][cur] = interval_length_vect[i][label][j];
	      interval_cpgCount[i][label][cur] =  interval_cpgCount[i][label][j];
	      interval_cpgMeth[i][label][cur] =  interval_cpgMeth[i][label][j];
	      cur++;
	    }
	    
	  }
	  interval_length_vect[i][label].resize(cur);
	  interval_cpgCount[i][label].resize(cur);
	  interval_cpgMeth[i][label].resize(cur);
	  /*
	  if (toDelete.size() > 0) {
	    didDelete=true;
	    cerr << "DELETED: " << toDelete.size() << endl;
	    cerr << "LENGTHS " << i << ":" ;
	    for (auto l: interval_length_vect[i][label]) {
	      cerr << " " << l;
	    }
	    cerr << endl;
	  }
	  */
	  //	  if (interval_length_vect[i][label][nValue-1] > loMean + 3*loSD + 100) {	  
	  //	    cout << label << " SUSPECT high length " << loMean << "\t" << loSD << "\t" << loMean + 3*loSD << "\t" << interval_length_vect[i][label][nValue-1] << endl;
	  //	  }
	}	  
      }

      for (auto h: {0,1}) {
	int totMeth=0, totCpg=0, totLen=0;
	for (int r=0; r < interval_cpgCount[h][label].size(); r++) {
	  totCpg += interval_cpgCount[h][label][r];
	  totMeth += interval_cpgMeth[h][label][r];
	  totLen += interval_length_vect[h][label][r];
	}
	float newAvgCpg, newAvgMeth;
	assert(interval_cpgCount[h][label].size() == interval_cpgMeth[h][label].size());
	if (interval_cpgCount[h][label].size() > 0) 
	  avgCpGCount[h] = ((float)totCpg) / interval_cpgCount[h][label].size();
	else
	  avgCpGCount[h] = -1;
	if (totCpg > 0)
	  avgMethRatio[h] = ((float)totMeth) / totCpg;
	else
	  avgMethRatio[h] = -1;
	if (interval_length_vect[h][label].size() > 0) 
	  avgLen[h] = ((float)totLen)/interval_length_vect[h][label].size();
	else
	  avgLen[h] = -1;
	
	/*
	if (didDelete > 0) {
	  cerr << "Prev cpg:  " << avgCpGCount[h] << " new: " << totCpg << endl;
	  cerr << "Prev meth:  " << avgMethRatio[h] << " new: " << newAvgMeth << endl;
	}
	*/
      }
      
      std::cout << interval_bed[label].chrom << "\t"
		<< interval_bed[label].start << "\t" 
		<< interval_bed[label].end << "\t"
		<< nReads[0] << "\t" << avgLen[0] << "\t" << avgCpGCount[0] << "\t" << avgMethRatio[0] << "\t"
		<< nReads[1] << "\t" << avgLen[1] << "\t" << avgCpGCount[1] << "\t" << avgMethRatio[1] << endl;
    }

    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    sam_close(in);
    return 0;
}

