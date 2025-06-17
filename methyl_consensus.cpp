#include <htslib/sam.h>
#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cstring>
#include <cstdlib>
#include <tuple>
#include <math.h>
#include "spoa/spoa.hpp"


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

      mod_positions.push_back(strPos);	
	/*	if (strand == 0) {
	  mod_positions.push_back(strPos);
	}
	else {
	  mod_positions.push_back(seq.size() - strPos - 1);
	}
	*/
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
    if (!hp_data) {
      //      cout << "WARNING, input data does not have HP phase tag" << endl;
      return -1;
    } 

    char type = hp_data[0];
    if (type == 'i') {
        found = true;
        return bam_aux2i(hp_data);
    } else if (type == 'I') {
        found = true;
        return static_cast<int>(bam_aux2i(hp_data));
    } else if (type == 'C') {
      return static_cast<int>(bam_aux2i(hp_data));
    }else {
        // Unexpected type
        return -1;
    }
}

int GetHapIndices(int hap, vector<int> &hapIdx) {
  hapIdx.clear();
  if (hap == 1) {
    hapIdx.push_back(0);
    return 0;
  }
  else if (hap == 2) {
    hapIdx.push_back(1);
    return 1;
  }
  else {
    hapIdx.push_back(0);
    hapIdx.push_back(1);
    return 2;
  }
}

struct AlignmentResult {
    string alignedA;
    string alignedB;
    vector<int> mapAtoB;  // For each position in A, its corresponding B index, or -1
    vector<int> mapBtoA;  // For each position in B, its corresponding A index, or -1
};

AlignmentResult needlemanWunsch(const string& A, const string& B, const string aMe="", const string bMe="") {
    int match = 1, mismatch = -1, gap = -1;
    int m = A.size(), n = B.size();

    vector<vector<int>> dp(m+1, vector<int>(n+1));
    vector<vector<char>> traceback(m+1, vector<char>(n+1));
    
    // Initialize DP table
    for (int i = 0; i <= m; ++i) {
        dp[i][0] = 0;
        traceback[i][0] = 'U';
    }
    for (int j = 0; j <= n; ++j) {
        dp[0][j] = 0;
        traceback[0][j] = 'L';
    }
    traceback[0][0] = '0';

    // Fill DP table
    if (aMe.size() == 0) {
      for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
	  int diag = dp[i-1][j-1] + (A[i-1] == B[j-1] ? match : mismatch);
	  int up = dp[i-1][j] + gap;
	  int left = dp[i][j-1] + gap;
	  dp[i][j] = max({diag, up, left});
	  if (dp[i][j] == diag) traceback[i][j] = 'D';
	  else if (dp[i][j] == up) traceback[i][j] = 'U';
	  else traceback[i][j] = 'L';
        }
      }
    }
    else {
      for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
	  int diag;
	  if (A[i-1] == B[j-1]) {
	    if (bMe[j-1] == '*') {
	      if (aMe[i-1] == '.') {
		diag = dp[i-1][j-1] + match + 2;
	      }
	      else if (aMe[i-1] == '*') {
		diag = dp[i-1][j-1] + match + 8;
	      }
	      else {
		diag = dp[i-1][j-1] + match;
	      }
	    }
	    else {
	      diag = dp[i-1][j-1] + match;
	    }
	  }
	  else {
	    diag = dp[i-1][j-1] + mismatch;
	  }
	  int up = dp[i-1][j] + gap;
	  int left = dp[i][j-1] + gap;
	  dp[i][j] = max({diag, up, left});
	  if (dp[i][j] == diag) traceback[i][j] = 'D';
	  else if (dp[i][j] == up) traceback[i][j] = 'U';
	  else traceback[i][j] = 'L';
        }
      }
    }

    // Traceback
    int i = m, j = n;
    string alignedA = "", alignedB = "";

    int maxScore=0;
    int maxi, maxj;
    for (int i=1; i <= m; i++) {
      if (dp[i][n] > maxScore) {
	maxScore=dp[i][n];
	maxi=i; maxj=n;
      }
    }
    for (int j=1; j <= n; j++) {
      if (dp[m][j] > maxScore) {
	maxScore=dp[m][j];
	maxi=m; maxj=j;
      }
    }
    
    vector<int> mapAtoB(m, -1);
    vector<int> mapBtoA(n, -1);

    int posA = maxi - 1, posB = maxj - 1;
    while (i > 0 || j > 0) {
        if (traceback[i][j] == 'D') {
            alignedA += A[i-1];
            alignedB += B[j-1];
            mapAtoB[i-1] = j-1;
            mapBtoA[j-1] = i-1;
            --i; --j;
        } else if (traceback[i][j] == 'U') {
            alignedA += A[i-1];
            alignedB += '-';
            --i;
        } else if (traceback[i][j] == 'L') {
            alignedA += '-';
            alignedB += B[j-1];
            --j;
        }
    }

    reverse(alignedA.begin(), alignedA.end());
    reverse(alignedB.begin(), alignedB.end());

    return {alignedA, alignedB, mapAtoB, mapBtoA};
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

void RevComp(string &seq, string&seq_rc) {	    
  seq_rc= seq;
  int n=seq.size();
  for (auto i = 0; i < seq.size(); i++) {
    if (seq[i] == 'A') { seq_rc[n-i-1] = 'T';}
    if (seq[i] == 'T') { seq_rc[n-i-1] = 'A';}
    if (seq[i] == 'C') { seq_rc[n-i-1] = 'G';}
    if (seq[i] == 'G') { seq_rc[n-i-1] = 'C';}	    
  }
}


int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: countme <input.bam> <input.bed> <sample_name>\n";
        return 1;
    }
    vector<string> labels;
    auto bed_intervals = read_bed_file(argv[2], labels);

    // map: interval label -> (sum of per-read avg methylation, count of reads contributing)
    vector<std::unordered_map<std::string, vector<string> >>  interval_seqs(2), interval_annot(2);    
    vector<std::unordered_map<std::string, vector< int> >>  interval_strands(2);    
    vector<std::unordered_map<std::string, vector<vector< int> > >>  interval_cpgs(2);
    vector<std::unordered_map<std::string, vector<vector< int> > >>  interval_methyl(2);
    vector<std::unordered_map<std::string, std::pair<int, int>> > interval_sums(2);
    vector<std::unordered_map<std::string, std::pair<int, int>> > interval_lengths(2);
    vector<std::unordered_map<std::string, vector<int> > > interval_length_vect(2);
    vector<std::unordered_map<std::string, vector<int> > > interval_cpgCount(2);
    vector<std::unordered_map<std::string, vector<int> > > interval_cpgMeth(2);    
    std::unordered_map<std::string, BedInterval> interval_bed;
    string bam_file = argv[1];
    string sample_name = argv[3];
    samFile* in = sam_open(argv[1], "r");
    if (!in) {
        std::cerr << "Error opening BAM file\n";
        return 1;
    }

    bam_hdr_t* hdr = sam_hdr_read(in);
    bam1_t* b = bam_init1();
    int proc_read=0;

    cerr << "Running mecons v0.1" << endl;
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
	  interval_seqs[h][bed.label()]=vector<string>();
	  interval_annot[h][bed.label()]=vector<string>();	  
	  interval_cpgs[h][bed.label()]=vector<vector<int> >();
	  interval_methyl[h][bed.label()]=vector<vector<int> >();
	  interval_strands[h][bed.label()]=vector<int>();	  
	}
      }
    }
    hts_idx_t* idx = sam_index_load(in, bam_file.c_str());

    for (auto& chrom: bed_intervals) {
      for (auto& bed: chrom.second) {
        int tid = bam_name2id(hdr, bed.chrom.c_str());

        hts_itr_t* iter = sam_itr_queryi(idx, tid, bed.start, bed.end);
        if (!iter) continue;
        while (sam_itr_next(in, iter, b) >= 0) {	

      
	  if (b->core.flag & BAM_FUNMAP) continue;
	  if (b->core.flag & BAM_FSUPPLEMENTARY) continue;		  
	  ++proc_read;
	  //	  if (proc_read == 10) {exit(0);}
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
	  int hi=GetHapIndices(haplotype, hapIndices);
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
	  string methString=seq;
	  string seq_rc;
	  if (strand == 1 ) {
	    RevComp(seq, seq_rc);
	    methString = seq_rc;
	  }

	  auto pos_c = parse_mm_tag(mm_str, methString, base, strand);

	  std::vector<int> read_to_ref(read_len, -1);
	  assert(ref_end - ref_start + 1 > 0);
	  std::vector<int> ref_to_read(ref_end - ref_start + 1, -1);	  
	  uint32_t* cigar = bam_get_cigar(b);
	  int read_pos = 0, ref_pos = b->core.pos;

	  for (uint32_t i = 0; i < b->core.n_cigar; ++i) {
            uint32_t op = bam_cigar_op(cigar[i]);
            uint32_t len = bam_cigar_oplen(cigar[i]);

            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
	      for (uint32_t j = 0; j < len; ++j) {
		if (strand == 0) {
		  if (read_pos < read_len) read_to_ref[read_pos++] = ref_pos++;
		}
		else {
		  assert(ref_pos - ref_start < ref_to_read.size());
		  ++read_pos;
		  if (read_pos > 0) read_to_ref[seq.size() - read_pos] = ref_pos++;
		}
		if (strand == 0) {
		  ref_to_read[ref_pos - ref_start] = read_pos;
		}
		else {
		  ref_to_read[ref_pos - ref_start] = seq.size() - read_pos - 1;
		}
	      }
            } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
	      for (uint32_t j=0; j< len; j++) {
		assert(ref_pos - ref_start < ref_to_read.size());
		  if (strand == 0) {
		    ref_to_read[ref_pos - ref_start] = read_pos;
		  }
		  else {
		    ref_to_read[ref_pos - ref_start] = seq.size() - read_pos - 1;
		  }

		ref_pos++;
	      }
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
	      readIntvEnd   = lower_bound(readAlnStart, readAlnEnd, bed.end+10);
	    }
	    else {
	      readIntvStart = lower_bound(readAlnStart, readAlnEnd, bed.start-10, NegComp);
	      readIntvEnd   = lower_bound(readAlnStart, readAlnEnd, bed.end+10, NegComp);
	    }

	    int refToReadIntvStart=-1;

	    if (bed.start >= ref_start and bed.start < ref_end) {
	      int idx=bed.start - ref_start;
	      refToReadIntvStart = ref_to_read[idx];
	      while (refToReadIntvStart == -1 and idx + 1 < ref_to_read.size()) {
		idx++;
		refToReadIntvStart = ref_to_read[idx];
	      }
	    }

	    int refToReadIntvEnd=-1;

	    if (bed.end >= ref_start and bed.end < ref_end) {
	      int idx=bed.end - ref_start;
	      refToReadIntvEnd = ref_to_read[idx];
	      while (refToReadIntvEnd == -1 and idx + 1 < ref_to_read.size()) {
		idx++;
		refToReadIntvEnd = ref_to_read[idx];
	      }
	    }

	    if ((strand == 0 and *readIntvStart > bed.start+10) or (strand == 1 and *readIntvEnd > bed.start+10)) {
	      //	      cout << "Read starting inside interval" << endl;
	      continue;
	    }
	    int readIntvStartIdx = readIntvStart - read_to_ref.begin()+10;
	    int readIntvEndIdx = readIntvEnd - read_to_ref.begin()-10;
	    //	    cout << "Intv match coordinates: " << strand << "\t" << readIntvStartIdx << "\t" << refToReadIntvStart << "\t" << readIntvEndIdx << "\t" << refToReadIntvEnd << endl;
	    
	    if (readIntvStart == read_to_ref.end() or readIntvEnd == read_to_ref.end() ) { continue;}
	    
	    if (readIntvStartIdx >= readIntvEndIdx) {
	    }

	    else {
	      nMatched++;
	    }
	    string trSeq;
	    int origReadIntvStartIdx=readIntvStartIdx;
	    int origReadIntvEndIdx=readIntvEndIdx;
	    trSeq=seq.substr(readIntvStartIdx, readIntvEndIdx - readIntvStartIdx);
	    /*
	    if (strand == 1) {
	      string revComp;
	      RevComp(trSeq, revComp);
	      trSeq = revComp;
	      }*/
	    vector<int> cpgPos, methPos;
	    string readSeq, readAnnot;
	    readSeq = trSeq;
            for (int i = readIntvStartIdx; i < readIntvEndIdx; ++i) {
	      int methIndex;
	      if (strand == 0) {
		methIndex =i;
	      }
	      else {
		methIndex =methString.size() - i - 2;
		if (methIndex < 0) {
		  methIndex = 0;
		}
	      }
	      if (( seq[i] == 'C' and seq[i+1] == 'G') || meth_status[methIndex] >= 127) {
		total_c++;
		cpgPos.push_back(i - readIntvStartIdx);
		if (meth_status[methIndex] >= 127) {
		  methPos.push_back(i - readIntvStartIdx);		      
		  ++methylated;
		  readAnnot.push_back('*');
		}
		else {
		  readAnnot.push_back('|');
		}
	      }
	      else {
		readAnnot.push_back(' ');
	      }
            }
	    /*
	    cout << "strand: " << (int) strand << endl;
	    cout << "hap:   " << hi << " " << trSeq << endl;
	    cout << "annot:   " << readAnnot << endl << endl;
	    */
	    //	    cout << "rcguess  " << seq.substr(seq.size() - readIntvEndIdx, readIntvEndIdx - readIntvStartIdx) << endl;
	    //	    cout << seq << endl;
	    for (auto h: hapIndices) {
	      interval_strands[h][bed.label()].push_back(strand);
	      interval_sums[h][bed.label()].first += methylated;
	      interval_sums[h][bed.label()].second += total_c;
	      //	      cout << read_name<< " hap: " << h << " adding length " << readIntvEndIdx - readIntvStartIdx<< " " << readIntvEndIdx  << " " << readIntvStartIdx << endl;
	      interval_seqs[h][bed.label()].push_back(trSeq);
	      interval_cpgs[h][bed.label()].push_back(cpgPos);
	      interval_annot[h][bed.label()].push_back(readAnnot);
	      interval_methyl[h][bed.label()].push_back(methPos);	    
	      interval_lengths[h][bed.label()].first += readIntvEndIdx - readIntvStartIdx;
	      interval_cpgCount[h][bed.label()].push_back(total_c);
	      interval_cpgMeth[h][bed.label()].push_back(methylated);
	      interval_length_vect[h][bed.label()].push_back(readIntvEndIdx - readIntvStartIdx);
	      interval_lengths[h][bed.label()].second += 1;
	    }
	    interval_bed[bed.label()] = bed;
	  }
	}
	//	cout << "For " << read_name << " ovp " << overlapping.size() << " matched: " << nMatched << endl;
      }
    }

    for (auto label: labels) {
      vector<double> avgMethRatio(2,0), avgCpGCount(2,0), avgLen(2,0), nReads(2,0);
      bool didDelete=false;
      for (auto i =0; i < 2; i++) {
	auto alignment_engine = spoa::AlignmentEngine::Create(
							      spoa::AlignmentType::kSW,  5, -4, -8, -6, -8, -6);

	spoa::Graph graph{};
	
	for (const auto& it : interval_seqs[i][label]) {
	  auto alignment = alignment_engine->Align(it, graph);
	  graph.AddAlignment(alignment, it);
	}
	string trConsensus;
	vector<uint32_t> consCov;
	trConsensus = graph.GenerateConsensus(interval_seqs[i][label].size()/3, &consCov);

	map<int,int> cpgCount, methCount, cpgCountRev, methCountRev;
	for (auto j = 1; j< interval_seqs[i][label].size(); j++ ) {
	  AlignmentResult alnRes = needlemanWunsch(trConsensus, interval_seqs[i][label][j]);
	  int bp=0;
	  string alnCpg="";
	  for ( auto ap=0; ap < alnRes.alignedB.size(); ap++) {
	    char mc=' ';
	    if (alnRes.alignedB[ap] != '-') {
	      for (auto cp: interval_cpgs[i][label][j]) {
		if (cp == bp) {
		  mc='|';
		  for (auto me: interval_methyl[i][label][j]) {
		    if (me ==  bp) {
		      mc = '*';
		    }
		  }
		  break;
		}
	      }
	      bp++;
	    }
	    else {
	      mc='-';
	    }
	    alnCpg.push_back(mc);
	  }


	  for (auto cpgPos: interval_cpgs[i][label][j]) {
	    int cpgDest = alnRes.mapBtoA[cpgPos];
	    if (cpgDest != -1) {
	      if (cpgCount.find(cpgDest) == cpgCount.end()) {
		cpgCount[cpgDest] = 0;
	      }
	      cpgCount[cpgDest]++;
	      if (methCount.find(cpgDest) == methCount.end()) {
		methCount[cpgDest] = 0;
	      }
	    }
	  }
	  for (auto methPos: interval_methyl[i][label][j]) {
	    int methDest = alnRes.mapBtoA[methPos];
	    if (methDest != -1) {
	      if (methCount.find(methDest) == methCount.end()) {
		methCount[methDest] = 0;
	      }
	      methCount[methDest]++;
	    }
	    // This shouldn't happen
	    if (cpgCount.find(methDest) == cpgCount.end()) {
	      cpgCount[methDest] = 1;
	    }
	  }
	}
	string meAnnot;
	for (auto c: trConsensus) {
	  meAnnot.push_back(' ');
	}
	for (auto cpgIter: cpgCount) {
	  if (cpgIter.second > 3) {
	    meAnnot[cpgIter.first] = '|';
	  }
	  float meFrac = ((float)methCount[cpgIter.first]) / cpgIter.second;
	  if (cpgIter.second >3 and meFrac > 0.25) {
	    meAnnot[cpgIter.first] = '.';
	  }
	  if (cpgIter.second >3 and  meFrac > 0.5) {
	    meAnnot[cpgIter.first] = '*';
	  }
	}
	/*
	  for (auto j = 1; j< interval_seqs[i][label].size(); j++ ) {
	  cout << "trAnnot: " << meAnnot << endl;
	  cout << "reAnnot: " << interval_annot[i][label][j] << endl << endl;
	  }
	*/


	for (auto j = 1; j< interval_seqs[i][label].size(); j++ ) {
	  AlignmentResult alnRes = needlemanWunsch(trConsensus, interval_seqs[i][label][j], meAnnot, interval_annot[i][label][j]);
	  int bp=0;
	  string alnCpg="";
	  for ( auto ap=0; ap < alnRes.alignedB.size(); ap++) {
	    char mc=' ';
	    if (alnRes.alignedB[ap] != '-') {
	      for (auto cp: interval_cpgs[i][label][j]) {
		if (cp == bp) {
		  mc='|';
		  for (auto me: interval_methyl[i][label][j]) {
		    if (me ==  bp) {
		      mc = '*';
		    }
		  }
		  break;
		}
	      }
	      bp++;
	    }
	    else {
	      mc='-';
	    }
	    alnCpg.push_back(mc);
	  }

	  for (auto cpgPos: interval_cpgs[i][label][j]) {
	    int cpgDest = alnRes.mapBtoA[cpgPos];
	    if (cpgDest != -1) {
	      if (cpgCountRev.find(cpgDest) == cpgCountRev.end()) {
		cpgCountRev[cpgDest] = 0;
	      }
	      cpgCountRev[cpgDest]++;
	      if (methCountRev.find(cpgDest) == methCountRev.end()) {
		methCountRev[cpgDest] = 0;
	      }
	    }
	  }
	  for (auto methPos: interval_methyl[i][label][j]) {
	    int methDest = alnRes.mapBtoA[methPos];
	    if (methDest != -1) {
	      if (methCountRev.find(methDest) == methCountRev.end()) {
		methCountRev[methDest] = 0;
	      }
	      methCountRev[methDest]++;
	    }
	    // This shouldn't happen
	    if (cpgCountRev.find(methDest) == cpgCountRev.end()) {
	      cpgCountRev[methDest] = 1;
	    }
	  }

	  
	}

	//	cout << "h" << i <<" " << trConsensus << endl;
	//	cout << "   " << meAnnot << endl;
	/*
	cout << sample_name << "\t" << label << "\t" << trConsensus.size() << "\t" << i << "\t";
	for (auto cpgIter: cpgCount) {
	  cout << cpgIter.first << "," << cpgIter.second << "," << ((float)methCount[cpgIter.first]) / cpgIter.second << "/";
	}
	cout << "-1,-1,-1" << endl;
	*/
	cout <<sample_name << "\t" << label << "\t" << trConsensus.size() << "\t" << i << "\t";
	for (auto cpgIter: cpgCountRev) {
	  cout << cpgIter.first << "," << cpgIter.second << "," << ((float)methCountRev[cpgIter.first]) / cpgIter.second << "/";
	}
	cout << "-1,-1,-1" << endl;
	
      	
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
	
      }
    }

    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    sam_close(in);
    return 0;
}

