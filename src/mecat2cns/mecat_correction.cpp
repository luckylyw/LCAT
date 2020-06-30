#include "mecat_correction.h"

#include "MECAT_AlnGraphBoost.H"

using namespace ns_banded_sw;

namespace ns_meap_cns {
//这里还要修改
#define FMAT 1
#define FDEL 2
#define FINS 4
#define UNDS 8 
#define M_FMAT 16
#define M_FDEL 32
#define M_FINS 64
#define M_UNDS 128

inline uint1
identify_one_consensus_item(CnsTableItem& cns_item, const int min_cov)
{
	uint1 ident = 0;
	int cov = cns_item.mat_cnt + cns_item.ins_cnt;
	if (cns_item.mat_cnt >= cov * 0.8) ident |= FMAT;
	if (cns_item.ins_cnt >= cov * 0.8) ident |= FINS;
	if (!ident) ident |= UNDS;
	if (cns_item.del_cnt >= cov * 0.4) ident |= FDEL;
	return ident;
}

inline uint1
identify_one_consensus_item(CnsTableItem& cns_item, CnsTableItem& cns_item_new, const int min_cov, const double min_modify_coverage)
{//修改了一下，存疑,还要修改一下阈值？？？？
	//cns_id_vec[i] = identify_one_consensus_item(cns_list[i], cns_list_new[i], min_cov);
	uint1 ident = 0;
	//std::cout << "ident : " << ident << "\n";
	//printf("printf<<ident:%c\n", ident);
	//uint1 struct_modified = 0;
	//std::cout << "mat_cnt=" << cns_item.mat_cnt <<";ins_cnt="<<cns_item.ins_cnt<<";del_cnt="<<cns_item.del_cnt<<";skip_cnt="<<cns_item.skip_cnt<<"\n";
	int cov_all = cns_item.mat_cnt + cns_item.ins_cnt;
	//std::cout << "cov_all=" << cov_all <<"&";
	int cov_part = cns_item_new.mat_cnt + cns_item_new.ins_cnt;
	//std::cout << "cov_part=" << cov_part <<"\n";

	//****加入的参数m:modify coverage threshold
	if (cns_item_new.skip_cnt >= cov_all * min_modify_coverage){//存在identity过低区域
		//std::cout << "***********modifiy    " ;
		if (cns_item_new.mat_cnt >= cov_part * 0.8) {
			ident |= M_FMAT;
			//std::cout << "M_FMAT&";
			//std::cout << "flag=" << !(ident&M_FMAT)<<"\n";
		}
		if (cns_item_new.ins_cnt >= cov_part * 0.8) {
			ident |= M_FINS;
			//std::cout << "M_FINS&";
			//std::cout << "flag=" << (ident&M_FINS)<<"\n";
		}
		if (!ident) {
			ident |= M_UNDS;
			//std::cout << "M_UNDS&";
			//std::cout << "flag=" << (ident&M_UNDS)<<"\n";
		}
		if (cns_item_new.del_cnt >= cov_part * 0.4) {
			ident |= M_FDEL;
			//std::cout << "M_FDEL&";
			//std::cout << "flag=" << (ident&M_FDEL)<<"\n";
		}
	}else{//不存在identity过低区域
		//std::cout << "===========normal    ";
		if (cns_item.mat_cnt >= cov_all * 0.8) {
			ident |= FMAT;
			//std::cout << "FMAT&";
			//std::cout << "flag=" << !(ident&FMAT)<<"\n";
		}
		if (cns_item.ins_cnt >= cov_all * 0.8) {
			ident |= FINS;
			//std::cout << "FINS&";
			//std::cout << "flag=" << (ident&FINS)<<"\n";
		}
		if (!ident) {
			ident |= UNDS;
			//std::cout << "UNDS&";
			//std::cout << "flag=" << (ident&UNDS)<<"\n";
		}
		if (cns_item.del_cnt >= cov_all * 0.4) {
			ident |= FDEL;
			//std::cout << "FDEL&";
			//std::cout << "flag=" << (ident&FDEL)<<"\n";
		}
	}

	return ident;
}

struct CompareOverlapByOverlapSize
{
	bool operator()(const Overlap& a ,const Overlap& b)
	{
		const index_t ovlp_a = std::max(a.qend - a.qoff, a.send - a.soff);
		const index_t ovlp_b = std::max(b.qend - b.qoff, b.send - b.soff);
		return ovlp_a > ovlp_b;
	}
};

void
meap_add_one_aln(const std::string& qaln, const std::string& saln, index_t start_soff, CnsTableItem* cns_table, const char* org_seq)
{//meap_add_one_aln(nqstr, ntstr, m5soff(*m5), cns_table, tstr.data());
	//meap_add_one_aln(nqstr, ntstr, m5soff(*m5), cns_table, tstr.data());
	r_assert(qaln.size() == saln.size());
	const index_t aln_size = qaln.size();//比对上的size
	index_t i = 0;
	const char kGap = '-';
	while (i < aln_size)
	{
		const char q = qaln[i];
		const char s = saln[i];
		if (q == kGap && s == kGap) { ++i; continue; }
		//修改的地方，存疑
		if(q == 'N' && s != '-'){
			++cns_table[start_soff].skip_cnt;
			cns_table[start_soff].base = s;
			++start_soff;
			++i;
		}else if( q == 'N' && s == '-'){
			//pass,不算delete错误
			index_t j = i + 1;
			while (j < aln_size && saln[j] == kGap) ++j;
			++cns_table[start_soff - 1].del_cnt;
			i = j;

		}else if(q == s){ ++cns_table[start_soff].mat_cnt; cns_table[start_soff].base = s; ++start_soff; ++i; }
		else if (q == kGap) { ++cns_table[start_soff].ins_cnt; ++start_soff; ++i; }
		else
		{
			r_assert(s == kGap);
			index_t j = i + 1;
			while (j < aln_size && saln[j] == kGap) ++j;
			++cns_table[start_soff - 1].del_cnt;
			i = j;
		}
	}
}
/*
//修改后
void
meap_add_one_aln(const std::string qaln, const std::string saln, index_t start_soff, CnsTableItem* cns_table, const char* org_seq)
{//meap_add_one_aln(nqstr, ntstr, m5soff(*m5), cns_table, tstr.data());
	//meap_add_one_aln(nqstr, ntstr, m5soff(*m5), cns_table, tstr.data());
	r_assert(qaln.size() == saln.size());
	const index_t aln_size = qaln.size();//比对上的size
	index_t i = 0;
	const char kGap = '-';
	while (i < aln_size)
	{
		const char q = qaln[i];
		const char s = saln[i];
		if (q == kGap && s == kGap) { ++i; continue; }
		
		if (q == s) { ++cns_table[start_soff].mat_cnt; cns_table[start_soff].base = s; ++start_soff; ++i; }
		else if (q == kGap) { ++cns_table[start_soff].ins_cnt; ++start_soff; ++i; }
		else
		{
			r_assert(s == kGap);
			index_t j = i + 1;
			while (j < aln_size && saln[j] == kGap) ++j;
			++cns_table[start_soff - 1].del_cnt;
			i = j;
		}
	}
}
*/

void
meap_cns_one_indel(const int sb, const int se, CnsAlns& cns_vec, 
				   const int min_cov, std::string& aux_qstr,
				   std::string& aux_tstr, std::string& cns)
{//meap_cns_one_indel(i + start_soff, j + start_soff, cns_vec, cns_list[i].mat_cnt + cns_list[i].ins_cnt, aux_qstr, aux_tstr, cns);
	int temp_min_cov = min_cov;//修改的地方
	AlnGraphBoost ag(se - sb + 1);
	int sb_out;//起始位置
	for (CnsAln* iter = cns_vec.begin(); iter != cns_vec.end(); ++iter)
	{
		if ((*iter).retrieve_aln_subseqs(sb, se, aux_qstr, aux_tstr, sb_out))//aux_qstr, qux_tstr是复杂区域片段，从cns_vec中复制得来的
		{
			ag.addAln(aux_qstr, aux_tstr, sb_out - sb + 1);
		}else{
			temp_min_cov--;
		}
	}		
	ag.mergeNodes();
	ag.consensus(temp_min_cov * 0.4, cns);//这个min_cov还需要再修改
}
void
meap_consensus_one_segment(CnsTableItem* cns_list, const int cns_list_size, 
						   uint1* cns_id_vec,
						   int start_soff, CnsAlns& cns_vec, 
						   std::string& aux_qstr, std::string& aux_tstr,
						   std::string& target, const int min_cov)
{
	for (int i = 0; i < cns_list_size; ++i) cns_id_vec[i] = identify_one_consensus_item(cns_list[i], min_cov);
	int i = 0, j; 
	std::string cns;
	target.clear();
	while (i < cns_list_size && !(cns_id_vec[i] & FMAT)) ++i;//非match时， ++i
	while (i < cns_list_size)
	{
		target.push_back(cns_list[i].base);
		j = i + 1;
		while (j < cns_list_size && !(cns_id_vec[j] & FMAT)) ++j;//非match时， ++j
		
		bool need_refinement = false;
		for (int k = i; k < j; ++k)
			if ((cns_id_vec[k] & UNDS) || (cns_id_vec[k] & FDEL)) { need_refinement = true; break; }
		if (need_refinement)
		{
			meap_cns_one_indel(i + start_soff, j + start_soff, cns_vec, cns_list[i].mat_cnt + cns_list[i].ins_cnt, aux_qstr, aux_tstr, cns);
			if (cns.size() > 2) target.append(cns.data() + 1, cns.size() - 2);
		}
		i = j;
	}
}

void
meap_consensus_one_segment(CnsTableItem* cns_list, CnsTableItem* cns_list_new,  const int cns_list_size, 
						   uint1* cns_id_vec,
						   int start_soff, CnsAlns& cns_vec, CnsAlns& cns_vec_new,
						   std::string& aux_qstr, std::string& aux_tstr,
						   std::string& target, const int min_cov, const double min_modify_coverage)
{//meap_consensus_one_segment(cns_table + beg, cns_table_new + beg, end - beg, id_list,
						//				   beg, cns_vec, cns_vec_new, aux_qstr, aux_tstr, cns_seq, min_cov);
	for (int i = 0; i < cns_list_size; ++i) {
		cns_id_vec[i] = identify_one_consensus_item(cns_list[i], cns_list_new[i], min_cov, min_modify_coverage);//判断这个地方的碱基是（MAT/INS/DEL/UNDEFINE）
		//std::cout << "cns_id_vec[" << i << "]=" << cns_id_vec[i] << " & ";
	}
	//std::cout << "\n";
	int i = 0, j; 
	std::string cns;
	target.clear();
	//改了这个地方的&&,找了好久才找到的错误，我太菜了，哭哭
	while (i < cns_list_size && (!(cns_id_vec[i] & FMAT) && !(cns_id_vec[i] & M_FMAT))) ++i;//如果不是match，++i
	//std::cout << "i" << i << "\n";
	//没有进入这个while, i=2325,实际上i应该=0
	while (i < cns_list_size)
	{
		//std::cout << "cns_list[" << i << "].base=" << cns_list[i].base << "\n";
		target.push_back(cns_list[i].base);//如果match，存放cns_table中这个位点的base
		//std::cout << target << "\n";
		j = i + 1;
		//同一个错误，
		while (j < cns_list_size && (!(cns_id_vec[j] & FMAT) && !(cns_id_vec[j] & M_FMAT))) ++j;//看后面的位点是不是match，不是，++j
		
		bool need_refinement_normal = false;
		bool need_refinement_modify = false;
		for (int k = i; k < j; ++k){
			if((cns_id_vec[k] & UNDS) || (cns_id_vec[k] & FDEL)){
				need_refinement_normal = true;
			}
			if((cns_id_vec[k] & M_UNDS) || (cns_id_vec[k] & M_FDEL)){
				need_refinement_modify = true;
				break;
			}
		}
		// 结构改变的地方，N算数
		if (need_refinement_modify)
		{//要修改的地方
			//修改，加入参数cns_vec_new
			//std::cout<<"need_refinement_modify"<<"\n";
			meap_cns_one_indel(i + start_soff, j + start_soff, cns_vec_new, cns_list_new[i].mat_cnt + cns_list_new[i].ins_cnt, aux_qstr, aux_tstr, cns);
			if (cns.size() > 2) target.append(cns.data() + 1, cns.size() - 2);
		}
		// 结构没有改变的地方，N不算数
		else if(need_refinement_normal)
		{
			//std::cout<<"need_refinement_normal"<<"\n";
			meap_cns_one_indel(i + start_soff, j + start_soff, cns_vec, cns_list[i].mat_cnt + cns_list[i].ins_cnt, aux_qstr, aux_tstr, cns);
			if (cns.size() > 2) target.append(cns.data() + 1, cns.size() - 2);
		}
		i = j;
	}
	//std::cout << "target:" << target << "\n";
}

struct CmpMappingRangeBySoff
{
	bool operator()(const MappingRange& a, const MappingRange& b)
	{
		return (a.start == b.start) ? (a.end > b.end) : (a.start < b.start);
	}
};

void 
get_effective_ranges(std::vector<MappingRange>& mranges, std::vector<MappingRange>& eranges, const int read_size, const int min_size)
{//get_effective_ranges(mranges, eranges, read_size, ctd->rco.min_size);
	eranges.clear();
	if (mranges.size() == 0) return;
	std::vector<MappingRange>::iterator iter;
	for (iter = mranges.begin(); iter != mranges.end(); ++iter)
		if (iter->start <= 500 && read_size - iter->end <= 500)
		{
			eranges.push_back(MappingRange(0, read_size));
			return;
		}

	std::sort(mranges.begin(), mranges.end(), CmpMappingRangeBySoff());
	const int nr = mranges.size();
	int i = 0, j;
	int left = mranges[i].start, right;
	while (i < nr)
	{
		j = i + 1;
		while (j < nr && mranges[j].end <= mranges[i].end) ++j;
		if (j == nr)
		{
			right = mranges[i].end;
			if (right - left >= min_size * 0.95) eranges.push_back(MappingRange(left, right));
			break;
		}
		if (mranges[i].end - mranges[j].start < 1000)
		{
			right = std::min(mranges[i].end, mranges[j].start);
			if (right - left >= min_size * 0.95) eranges.push_back(MappingRange(left, right));
			left = std::max(mranges[i].end, mranges[j].start);
		}
		i = j;
	}
}

void
output_cns_result(std::vector<CnsResult>& cns_results,
				  CnsResult& cr,
				  const index_t beg,
				  const index_t end,
				  std::string& cns_seq) 
{
	const size_t MaxSeqSize = 60000;
	const size_t OvlpSize = 10000;
	// BlkSize must be >= OvlpSize
	const size_t BlkSize = MaxSeqSize - OvlpSize - 1000;
	
	const size_t size = cns_seq.size();
	if (size <= MaxSeqSize) {
		cr.range[0] = beg;
		cr.range[1] = end;
		cr.seq = cns_seq;
		cns_results.push_back(cr);//id,range,se
		//std::cout<<"cr.seq:"<<cr.seq<<"\n";
	} else {
		const size_t cutoff = size - OvlpSize - 1000;
		size_t L = 0, R;
		do {
			R = L + BlkSize;
			if (R >= cutoff) {
				R = size;
			}
			cr.range[0] = L + beg;
			cr.range[1] = R < size && R + beg < static_cast<size_t>(end) ? R + beg : end;
			cr.seq = cns_seq.substr(L, R - L);
			cns_results.push_back(cr);
			L = R - OvlpSize;
			//std::cout<<"cr.seq:"<<cr.seq<<"\n";
		} while (R < size);
	}
}

inline bool
check_ovlp_mapping_range(const int qb, const int qe, const int qs,
						 const int sb, const int se, const int ss,
						 double ratio)
{
	const int oq = qe - qb;
	const int qqs = qs * ratio;
	const int os = se - sb;
	const int qss = ss * ratio;
	return oq >= qqs || os >= qss;
}
void
consensus_worker(CnsTableItem* cns_table,
				 uint1* id_list,
				 CnsAlns& cns_vec,
				 std::string& aux_qstr,
				 std::string& aux_tstr,
				 std::vector<MappingRange>& eranges,
				 const int min_cov,
				 const int min_size,
				 const int read_id,
				 std::vector<CnsResult>& cns_results)
{
	index_t beg = 0, end;
	CnsResult cns_result;
	std::string cns_seq;
	cns_result.id = read_id;
	std::vector<MappingRange>::iterator miter;
	for (miter = eranges.begin(); miter != eranges.end(); ++miter)
	{
		int L = miter->start, R = miter->end;
		beg = L;
		while (beg < R)
		{
			while (beg < R && cns_table[beg].mat_cnt + cns_table[beg].ins_cnt < min_cov) ++beg;
			end = beg + 1;
			while (end < R && cns_table[end].mat_cnt + cns_table[end].ins_cnt >= min_cov) ++end;
			if (end - beg >= 0.95 * min_size)
			{
				meap_consensus_one_segment(cns_table + beg, end - beg, id_list,
										   beg, cns_vec, aux_qstr, aux_tstr, cns_seq, min_cov);
				
				if (cns_seq.size() >= min_size) output_cns_result(cns_results, cns_result, beg, end, cns_seq);
			}
			
			beg = end;
		}
	}
}
//修改，加入了参数cns_vec_new, cns_table_new
void
consensus_worker(CnsTableItem* cns_table, CnsTableItem* cns_table_new,
				 uint1* id_list,
				 CnsAlns& cns_vec,
				 CnsAlns& cns_vec_new,
				 std::string& aux_qstr,
				 std::string& aux_tstr,
				 std::vector<MappingRange>& eranges,
				 const int min_cov,
				 const int min_size,
				 const int read_id,
				 std::vector<CnsResult>& cns_results,
				 const double min_modify_coverage)
{//consensus_worker(cns_table, cns_table_new, ctd->id_list, cns_vec, cns_vec_new, nqstr, ntstr, eranges, ctd->rco.min_cov, ctd->rco.min_size, read_id,  cns_results);
	//std::cout << "enter consensus_worker!!!" << "\n";
	index_t beg = 0, end;
	CnsResult cns_result;//id,reange[2],seq
	std::string cns_seq;
	cns_result.id = read_id;
	std::vector<MappingRange>::iterator miter;//start,end
	for (miter = eranges.begin(); miter != eranges.end(); ++miter)
	{
		int L = miter->start, R = miter->end;
		beg = L;
		while (beg < R)
		{
			//找到覆盖度(mat_cnt + ins_cnt)>min_cov的区域，起始base位点分别为beg,end
			//修改了这里，加入了cns_table[beg].skip_cnt
			while (beg < R && cns_table[beg].mat_cnt + cns_table[beg].ins_cnt + cns_table[beg].skip_cnt < min_cov) ++beg;
			end = beg + 1;
			while (end < R && cns_table[end].mat_cnt + cns_table[end].ins_cnt + cns_table[beg].skip_cnt >= min_cov) ++end;
			if (end - beg >= 0.95 * min_size)
			{
				//std::cout << "beg:" << beg << "#end:" << end << "\n";
				//修改，加入了参数cns_vec_new,加入参数cns_table_new
				meap_consensus_one_segment(cns_table + beg, cns_table_new + beg, end - beg, id_list,
										   beg, cns_vec, cns_vec_new, aux_qstr, aux_tstr, cns_seq, min_cov, min_modify_coverage);
				//std::cout << "cns_seq:" << cns_seq << "\n";
				//把生成的cns_seq赋值给cns_result，并push进入cns_results.
				if (cns_seq.size() >= min_size) output_cns_result(cns_results, cns_result, beg, end, cns_seq);
			}
			
			beg = end;
		}
	}
}

void
consensus_one_read_m4_pacbio(ConsensusThreadData* ctd, const index_t read_id, const index_t sid, const index_t eid)
{
	PackedDB& reads = *ctd->reads;
	ExtensionCandidate* overlaps = ctd->candidates;
	DiffRunningData* drd_s = ctd->drd_s;
	DiffRunningData* drd = NULL;
	M5Record* m5 = ctd->m5;
	CnsAlns& cns_vec = ctd->cns_alns;
	std::vector<CnsResult>& cns_results = ctd->cns_results;
	const index_t read_size = overlaps[sid].ssize;
	std::vector<char>& qstr = ctd->query;
	std::vector<char>& tstr = ctd->target;
	tstr.resize(read_size);
	reads.GetSequence(read_id, true, tstr.data(), read_size);
	std::string& nqstr = ctd->qaln;
	std::string& ntstr = ctd->saln;
	const int min_align_size = ctd->rco.min_align_size;
	const int max_added = 60;

	index_t L, R;
	if (eid - sid <= max_added)
	{
		L = sid;
		R = eid;
	}
	else
	{
		L = sid;
		R = L + max_added;    
		std::sort(overlaps + sid, overlaps + eid, CompareOverlapByOverlapSize());
	}

	CnsTableItem* cns_table = ctd->cns_table;
	std::for_each(cns_table, cns_table + read_size, CnsTableItemCleaner());
	cns_vec.clear();
	for (index_t i = L; i < R; ++i)
	{
		Overlap& ovlp = overlaps[i];
		qstr.resize(ovlp.qsize);
		reads.GetSequence(ovlp.qid, ovlp.qdir == FWD, qstr.data(), ovlp.qsize);
		index_t qext = ovlp.qext;
		index_t sext = ovlp.sext;
		if (ovlp.qdir == REV) qext = ovlp.qsize - 1 - qext;
		drd = drd_s;
		bool r = GetAlignment(qstr.data(), qext, qstr.size(), tstr.data(), sext, tstr.size(), drd, *m5, 0.15, min_align_size);
		if (r)
		{
			normalize_gaps(m5qaln(*m5), m5saln(*m5), strlen(m5qaln(*m5)), nqstr, ntstr, true);
			meap_add_one_aln(nqstr, ntstr, m5soff(*m5), cns_table, tstr.data());
			cns_vec.add_aln(m5soff(*m5), m5send(*m5), nqstr, ntstr);
		}
	}
	
	std::vector<MappingRange> mranges, eranges;
	cns_vec.get_mapping_ranges(mranges);
	get_effective_ranges(mranges, eranges, read_size, ctd->rco.min_size);

	consensus_worker(cns_table, ctd->id_list, cns_vec, nqstr, ntstr, eranges, ctd->rco.min_cov, ctd->rco.min_size, read_id,  cns_results);
}

void
consensus_one_read_m4_nanopore(ConsensusThreadData* ctd, const index_t read_id, const index_t sid, const index_t eid)
{
	PackedDB& reads = *ctd->reads;
	ExtensionCandidate* overlaps = ctd->candidates;
	DiffRunningData* drd_s = ctd->drd_s;
	DiffRunningData* drd = NULL;
	M5Record* m5 = ctd->m5;
	CnsAlns& cns_vec = ctd->cns_alns;
	std::vector<CnsResult>& cns_results = ctd->cns_results;
	const index_t read_size = overlaps[sid].ssize;
	std::vector<char>& qstr = ctd->query;
	std::vector<char>& tstr = ctd->target;
	tstr.resize(read_size);
	reads.GetSequence(read_id, true, tstr.data(), read_size);
	std::string& nqstr = ctd->qaln;
	std::string& ntstr = ctd->saln;
	const int min_align_size = ctd->rco.min_align_size;
	const double min_mapping_ratio = ctd->rco.min_mapping_ratio - 0.02;

	index_t L, R;
	if (eid - sid <= MAX_CNS_OVLPS)
	{
		L = sid;
		R = eid;
	}
	else
	{
		L = sid;
		R = L + MAX_CNS_OVLPS;    
		std::sort(overlaps + sid, overlaps + eid, CompareOverlapByOverlapSize());
	}

	CnsTableItem* cns_table = ctd->cns_table;
	std::for_each(cns_table, cns_table + read_size, CnsTableItemCleaner());
	cns_vec.clear();
	for (index_t i = L; i < R; ++i)
	{
		Overlap& ovlp = overlaps[i];
		qstr.resize(ovlp.qsize);
		reads.GetSequence(ovlp.qid, ovlp.qdir == FWD, qstr.data(), ovlp.qsize);
		index_t qext = ovlp.qext;
		index_t sext = ovlp.sext;
		if (ovlp.qdir == REV) qext = ovlp.qsize - 1 - qext;
		drd = drd_s;
		bool r = GetAlignment(qstr.data(), qext, qstr.size(), tstr.data(), sext, tstr.size(), drd, *m5, 0.20, min_align_size);
		if (r && check_ovlp_mapping_range(m5qoff(*m5), m5qend(*m5), ovlp.qsize, m5soff(*m5), m5send(*m5), ovlp.ssize, min_mapping_ratio))
		{
			normalize_gaps(m5qaln(*m5), m5saln(*m5), strlen(m5qaln(*m5)), nqstr, ntstr, true);
			meap_add_one_aln(nqstr, ntstr, m5soff(*m5), cns_table, tstr.data());
			cns_vec.add_aln(m5soff(*m5), m5send(*m5), nqstr, ntstr);
		}
	}
	
	std::vector<MappingRange> mranges, eranges;
	eranges.push_back(MappingRange(0, read_size));

	consensus_worker(cns_table, ctd->id_list, cns_vec, nqstr, ntstr, eranges, ctd->rco.min_cov, ctd->rco.min_size, read_id,  cns_results);
}

struct CmpExtensionCandidateByScore
{
	bool operator()(const ExtensionCandidate& a, const ExtensionCandidate& b)
	{
        if (a.score != b.score) return a.score > b.score;
		if (a.qid != b.qid) return a.qid < b.qid;
		return a.qext < b.qext;
	}
};

inline bool
check_cov_stats(u1_t* cov_stats, int soff, int send)
{
	const int max_cov = 20;
	int n = 0;
	for (int i = soff; i < send; ++i)
		if (cov_stats[i] >= max_cov)
			++n;
	if (send - soff >= n + 200) 
	{
		for (int i = soff; i < send; ++i) ++cov_stats[i];
		return true;
	}
	return false;
}

void
consensus_one_read_can_pacbio(ConsensusThreadData* ctd, const index_t read_id, const index_t sid, const index_t eid)
{
	PackedDB& reads = *ctd->reads;
	ExtensionCandidate* candidates = ctd->candidates;
	DiffRunningData* drd_s = ctd->drd_s;
	DiffRunningData* drd = NULL;
	M5Record* m5 = ctd->m5;
	CnsAlns& cns_vec = ctd->cns_alns;
	//std::cout<<"ctd->cns_alns="<<ctd->cns_alns<<std::endl;
	//std::cout<<"cns_vec="<<&cns_vec<<std::endl;
	//修改,加入的
	CnsAlns& cns_vec_new = ctd->new_cns_alns;//加入的
	std::vector<CnsResult>& cns_results = ctd->cns_results;
	const index_t read_size = candidates[sid].ssize;
	std::vector<char>& qstr = ctd->query;//候选序列
	std::vector<char>& tstr = ctd->target;//目标序列
	tstr.resize(read_size);
	reads.GetSequence(read_id, true, tstr.data(), read_size);
	std::string& nqstr = ctd->qaln;//候选比对结果
	std::string& ntstr = ctd->saln;//目标比对结果
	//后加入的
	std::string newnqstr;//加入的
	std::string& newnqstr2 = newnqstr;//加入的
	const int min_align_size = ctd->rco.min_align_size;
	const double min_mapping_ratio = ctd->rco.min_mapping_ratio - 0.02;
	const int max_added = 60;
	
	std::sort(candidates + sid, candidates + eid, CmpExtensionCandidateByScore());
	int num_added = 0;
    int num_ext = 0;
    const int max_ext = 200;
	CnsTableItem* cns_table = ctd->cns_table;
	CnsTableItem cns_table1[read_size];//加入的
	memcpy(cns_table1, cns_table, read_size);//加入的
	CnsTableItem* cns_table_new = cns_table1;//加入的
	std::for_each(cns_table_new, cns_table_new + read_size, CnsTableItemCleaner());//加入的
	std::for_each(cns_table, cns_table + read_size, CnsTableItemCleaner());
	cns_vec.clear();
	//修改，加入的
	cns_vec_new.clear();//加入的
	std::set<int> used_ids;
	u1_t* cov_stats = ctd->id_list;
	std::fill(cov_stats, cov_stats + MAX_SEQ_SIZE, 0);
	for (idx_t i = sid; i < eid && num_added < max_added && num_ext < max_ext; ++i)
	{
        ++num_ext;
		ExtensionCandidate& ec = candidates[i];
		r_assert(ec.sdir == FWD);
		if (used_ids.find(ec.qid) != used_ids.end()) continue;
		qstr.resize(ec.qsize);
		reads.GetSequence(ec.qid, ec.qdir == FWD, qstr.data(), ec.qsize);
		index_t qext = ec.qext;
		index_t sext = ec.sext;
		if (ec.qdir == REV) qext = ec.qsize - 1 - qext;
		drd = drd_s;
		bool r = GetAlignment(qstr.data(), qext, qstr.size(), tstr.data(), sext, tstr.size(), drd, *m5, 0.15, min_align_size);
		if (r && check_ovlp_mapping_range(m5qoff(*m5), m5qend(*m5), ec.qsize, m5soff(*m5), m5send(*m5), ec.ssize, min_mapping_ratio))
		{
			if (check_cov_stats(cov_stats, m5soff(*m5), m5send(*m5)))
			{
				++num_added;
				used_ids.insert(ec.qid);
				normalize_gaps(m5qaln(*m5), m5saln(*m5), strlen(m5qaln(*m5)), nqstr, ntstr, true);
				//加入的
				//****加入的参数:slide window length
    			const int k = ctd->rco.window_length;
    			//****加入的参数:identity threshold
    			const double identity_threshold = ctd->rco.identity_threshold;
				//int k = 75;
    			//double identity_threshold = 0.65;
    			slide_window2(nqstr, ntstr, newnqstr2, k, identity_threshold);

				meap_add_one_aln(nqstr, ntstr, m5soff(*m5), cns_table, tstr.data());
				meap_add_one_aln(newnqstr2, ntstr, m5soff(*m5), cns_table_new, tstr.data());

				cns_vec_new.add_aln(m5soff(*m5), m5send(*m5), newnqstr2, ntstr);
				cns_vec.add_aln(m5soff(*m5), m5send(*m5), nqstr, ntstr);
			}
		}
	}
	
	/*
	std::vector<MappingRange> mranges, eranges;
	std::cout << "std::vector<MappingRange> mranges, eranges; end" << "\n";
	eranges.push_back(MappingRange(0, read_size));
	std::cout << "eranges.push_back(MappingRange(0, read_size)); end" << "\n";
	consensus_worker(cns_table, cns_table_new, ctd->id_list, cns_vec, cns_vec_new, nqstr, ntstr, eranges, ctd->rco.min_cov, ctd->rco.min_size, read_id,  cns_results);
	*/
	std::vector<MappingRange> mranges, eranges;
	cns_vec.get_mapping_ranges(mranges);
	get_effective_ranges(mranges, eranges, read_size, ctd->rco.min_size);
	//****加入的参数:min_modify_coverage
    const double min_modify_coverage = ctd->rco.min_modify_coverage;
	consensus_worker(cns_table, cns_table_new, ctd->id_list, cns_vec, cns_vec_new, nqstr, ntstr, eranges, ctd->rco.min_cov, ctd->rco.min_size, read_id,  cns_results, min_modify_coverage);
}
/*
void
consensus_one_read_can_nanopore(ConsensusThreadData* ctd, const index_t read_id, const index_t sid, const index_t eid)
{
	//entry,ns_meap_cns::consensus_one_read_can_nanopore(&cns_data, sid, i, j)
	//entry，sid/readId为模板读数id，i/sid为模板在num_candidates上的位置，j/eid为候选在num_candidates上的位置？？
	//sid = candidates[i].sid;
	std::cout << "--------consensus_one_read_can_nanopore::start--------" << "\n";
	std::cout << "consensus_one_read_can_nanopore(ConsensusThreadData* ctd, const index_t read_id, const index_t sid, const index_t eid):::\n"<< "read_id = " << read_id << "\n";//模板读数id
	std::cout << "sid = " << sid << "\n";//第一个候选在num_candidates上的位置
	std::cout << "eid = " << eid << "\n";//最后一个候选在num_candidates上的位置
	PackedDB& reads = *ctd->reads;
	//std::cout << "ctd->reads = reads.id = " << reads->id << "\t reads.size = " << reads->size << "\n";
	ExtensionCandidate* candidates = ctd->candidates;
	DiffRunningData* drd_s = ctd->drd_s;
	DiffRunningData* drd = NULL;
	M5Record* m5 = ctd->m5;
	CnsAlns& cns_vec = ctd->cns_alns;
	std::vector<CnsResult>& cns_results = ctd->cns_results;
	const index_t read_size = candidates[sid].ssize;
	std::vector<char>& qstr = ctd->query;
	std::vector<char>& tstr = ctd->target;
	tstr.resize(read_size);
	reads.GetSequence(read_id, true, tstr.data(), read_size);
	std::string& nqstr= ctd->qaln;//引用
	std::string& ntstr = ctd->saln;//引用
	//std::cout << "nqstr=ctd->qaln"<<nqstr<<"\n";//无
	//std::cout << "ntstr=ctd->saln"<<ntstr<<"\n";//无
	const int min_align_size = ctd->rco.min_align_size;
	const double min_mapping_ratio = ctd->rco.min_mapping_ratio - 0.02;
	
	std::sort(candidates + sid, candidates + eid, CmpExtensionCandidateByScore());
	int num_added = 0;
    int num_ext = 0;
    const int max_ext = 200;
	CnsTableItem* cns_table = ctd->cns_table;
	std::for_each(cns_table, cns_table + read_size, CnsTableItemCleaner());
	cns_vec.clear();
	std::set<int> used_ids;
	u1_t* cov_stats = ctd->id_list;
	std::fill(cov_stats, cov_stats + MAX_SEQ_SIZE, 0);
	for (idx_t i = sid; i < eid && num_added < MAX_CNS_OVLPS && num_ext < max_ext; ++i)
	{//遍历每一个候选
		std::cout<<"start:::::i = " << i << "\n";
        ++num_ext;
		ExtensionCandidate& ec = candidates[i];
		std::cout << "candidates[i].qid = " << ec.qid << "\n";
		r_assert(ec.sdir == FWD);
		if (used_ids.find(ec.qid) != used_ids.end()) continue;
		qstr.resize(ec.qsize);
		reads.GetSequence(ec.qid, ec.qdir == FWD, qstr.data(), ec.qsize);
		index_t qext = ec.qext;
		index_t sext = ec.sext;
		if (ec.qdir == REV) qext = ec.qsize - 1 - qext;
		drd = drd_s;
		//std::cout<<"qstr.data() = "<<qstr.data()<<"\n";//ec.qid，无
		//std::cout<<"tstr.data() = "<<tstr.data()<<"\n";//模板读数数据，无
		std::cout<<"GetAlignment start"<<"\n";
		bool r = GetAlignment(qstr.data(), qext, qstr.size(), tstr.data(), sext, tstr.size(), drd, *m5, 0.20, min_align_size);//entry
		//std::cout<<"qstr.data() = "<<qstr.data()<<"\n";//无
		//std::cout<<"tstr.data() = "<<tstr.data()<<"\n";//无
		//std::cout << "m5qid(*m5) = " << m5qid(*m5) << "\t m5sid(*m5)" << m5sid(*m5) << "\n";
		if (r && check_ovlp_mapping_range(m5qoff(*m5), m5qend(*m5), ec.qsize, m5soff(*m5), m5send(*m5), ec.ssize, min_mapping_ratio))
		{
			if (check_cov_stats(cov_stats, m5soff(*m5), m5send(*m5)))
			{
				++num_added;
				used_ids.insert(ec.qid);
				std::cout << "--------normalize_gaps start.--------" << "\n";
				normalize_gaps(m5qaln(*m5), m5saln(*m5), strlen(m5qaln(*m5)), nqstr, ntstr, true);
				std::cout << "nqstr : " << nqstr << "\n";
    			std::cout << "ntstr : " << ntstr << "\n";
    			std::cout << "--------normalize_gaps end.--------" << "\n";
    			std::cout << "--------slide_window start.--------" << "\n";
    			int k = 30;
				slide_window(nqstr, ntstr, k);
				std::cout << "--------slide_window end.--------" << "\n";
				meap_add_one_aln(nqstr, ntstr, m5soff(*m5), cns_table, tstr.data());
				cns_vec.add_aln(m5soff(*m5), m5send(*m5), nqstr, ntstr);
			}
		}
		std::cout<<"end:::::i = " << i << "\n";
	}
	std::vector<MappingRange> mranges, eranges;
	eranges.push_back(MappingRange(0, read_size));
	consensus_worker(cns_table, ctd->id_list, cns_vec, nqstr, ntstr, eranges, ctd->rco.min_cov, ctd->rco.min_size, read_id,  cns_results);
	std::cout<<"--------consensus_one_read_can_nanopore::end--------"<<"\n";
}

}// namespace ns_meap_cns 

*/

void

consensus_one_read_can_nanopore(ConsensusThreadData* ctd, const index_t read_id, const index_t sid, const index_t eid)
//ns_meap_cns::consensus_one_read_can_nanopore(&cns_data, sid, i, j);//entry
{
	//std::cout << "--------consensus_one_read_can_nanopore::start--------" << "\n";
	//std::cout << "consensus_one_read_can_nanopore(ConsensusThreadData* ctd, const index_t read_id, const index_t sid, const index_t eid):::\n"<< "read_id = " << read_id << "\n";//模板读数id
	//std::cout << "sid = " << sid << "\n";//第一个候选在num_candidates上的位置
	//std::cout << "eid = " << eid << "\n";//最后一个候选在num_candidates上的位置
	PackedDB& reads = *ctd->reads;
	ExtensionCandidate* candidates = ctd->candidates;
	DiffRunningData* drd_s = ctd->drd_s;
	DiffRunningData* drd = NULL;
	M5Record* m5 = ctd->m5;
	CnsAlns& cns_vec = ctd->cns_alns;
	//修改,加入的
	CnsAlns& cns_vec_new = ctd->new_cns_alns;//加入的
	std::vector<CnsResult>& cns_results = ctd->cns_results;
	const index_t read_size = candidates[sid].ssize;
	std::vector<char>& qstr = ctd->query;
	std::vector<char>& tstr = ctd->target;
	tstr.resize(read_size);
	reads.GetSequence(read_id, true, tstr.data(), read_size);
	std::string& nqstr = ctd->qaln;
	std::string& ntstr = ctd->saln;
	//后加入的
	std::string newnqstr;//加入的
	std::string& newnqstr2 = newnqstr;//加入的
	const int min_align_size = ctd->rco.min_align_size;
	const double min_mapping_ratio = ctd->rco.min_mapping_ratio - 0.02;
	std::sort(candidates + sid, candidates + eid, CmpExtensionCandidateByScore());
	int num_added = 0;
    int num_ext = 0;
    const int max_ext = 200;
    //修改的，后面加入的
	CnsTableItem* cns_table = ctd->cns_table;
	CnsTableItem cns_table1[read_size];//加入的
	memcpy(cns_table1, cns_table, read_size);//加入的
	CnsTableItem* cns_table_new = cns_table1;//加入的
	std::for_each(cns_table_new, cns_table_new + read_size, CnsTableItemCleaner());//加入的

	std::for_each(cns_table, cns_table + read_size, CnsTableItemCleaner());
	cns_vec.clear();
	//修改，加入的
	cns_vec_new.clear();//加入的
	std::set<int> used_ids;
	u1_t* cov_stats = ctd->id_list;
	std::fill(cov_stats, cov_stats + MAX_SEQ_SIZE, 0);
	for (idx_t i = sid; i < eid && num_added < MAX_CNS_OVLPS && num_ext < max_ext; ++i)//这有个过滤，一点点，1.overlap条数不过多；2.候选条数不过多
	{
        ++num_ext;
		ExtensionCandidate& ec = candidates[i];
		/*
		std::string newnqstr;
		std::string& newnqstr2 = newnqstr;
		*/
		/*
qdir, qid, qext, qsize, qoff, qend;
	int sdir, sid, sext, ssize, soff, send;
		*/
		//std::cout << "candidates[i].qid = " << ec.qid << "\n";
		//std::cout << "candidates[i].qext = " << ec.qext << "\n";
		//std::cout << "candidates[i].qsize = " << ec.qsize << "\n";
		//std::cout << "candidates[i].qoff = " << ec.qoff << "\n";
		//std::cout << "candidates[i].qend = " << ec.qend << "\n";
		//std::cout << "candidates[i].sid = " << ec.sid << "\n";
		//std::cout << "candidates[i].sext = " << ec.sext << "\n";
		//std::cout << "candidates[i].ssize = " << ec.ssize << "\n";
		//std::cout << "candidates[i].soff = " << ec.soff << "\n";
		//std::cout << "candidates[i].send = " << ec.send << "\n";
		r_assert(ec.sdir == FWD);
		if (used_ids.find(ec.qid) != used_ids.end()) continue;//存疑，看不懂这句话啊
		qstr.resize(ec.qsize);
		reads.GetSequence(ec.qid, ec.qdir == FWD, qstr.data(), ec.qsize);
		index_t qext = ec.qext;
		index_t sext = ec.sext;
		if (ec.qdir == REV) qext = ec.qsize - 1 - qext;
		drd = drd_s;
		//std::cout<<"--------GetAlignment start--------"<<"\n";
		//std::cout << "qext = " << qext << "\n";
		//std::cout << "sext = " << sext << "\n";
		//std::cout << "qstr.size() = " << qstr.size() << "\n";
		//std::cout << "tstr.size() = " << tstr.size() << "\n";
		//std::cout << "qstr.data() = " << qstr.data() << "\n";
		//std::cout << "tstr.data() = " << tstr.data() << "\n";
		//std::cout << "GetAlignment()" << "\n";
		bool r = GetAlignment(qstr.data(), qext, qstr.size(), tstr.data(), sext, tstr.size(), drd, *m5, 0.20, min_align_size);
		/*
		std::cout << "qstr:";
		for (index_t i = 0; i < strlen(m5qaln(*m5)); ++i)
    	{
    		std::cout<<m5qaln(*m5)[i];
    	}
    	std::cout << "\n";
    	std::cout << "tstr:";
    	for (index_t i = 0; i < strlen(m5qaln(*m5)); ++i)
    	{
    		std::cout<<m5saln(*m5)[i];
    	}
    	std::cout << "\n";
		std::cout<<"--------GetAlignmetn end--------"<<"\n";
		*/
		if (r && check_ovlp_mapping_range(m5qoff(*m5), m5qend(*m5), ec.qsize, m5soff(*m5), m5send(*m5), ec.ssize, min_mapping_ratio))
		//if (r)
		{
			if (check_cov_stats(cov_stats, m5soff(*m5), m5send(*m5)))
			{
				++num_added;
				used_ids.insert(ec.qid);
				//std::cout << "--------normalize_gaps start.--------" << "\n";
				normalize_gaps(m5qaln(*m5), m5saln(*m5), strlen(m5qaln(*m5)), nqstr, ntstr, true);
				/*
				std::cout << "nqstr_id:" << ec.qid << "\n";
				std::cout << "nqstr_size=ec.qsize :" << ec.qsize << "\n";
				std::cout << "nqstr_size = m5qsize(*m5) :" << m5qsize(*m5) << "\n";
				std::cout << "m5qoff(m5)= drd->result->query_start + qrb = " << m5qoff(*m5) << "\n";
				std::cout << "m5qend(m5) = drd->result->query_end - qre = " << m5qend(*m5) << "\n";
				std::cout << "m5qdir(m5) = FWD = " << m5qdir(*m5) << "\n";
				std::cout << "nqstr : " << nqstr << "\n";

				std::cout << "ntstr_id:" << ec.sid << "\n";
				std::cout << "ntstr_size=ec.ssize :" << ec.ssize << "\n";
				std::cout << "ntstr_size = m5ssize(*m5) :" << m5ssize(*m5) << "\n";
				std::cout << "m5toff(m5)= drd->result->target_start + trb = " << m5soff(*m5) << "\n";
				std::cout << "m5tend(m5) = drd->result->target_end - tre = " << m5send(*m5) << "\n";
				std::cout << "m5tdir(m5) = FWD = " << m5sdir(*m5) << "\n";
    			std::cout << "ntstr : " << ntstr << "\n";
    			*/
    			//std::cout << ec.qid << " " << m5qoff(*m5) << " " << nqstr << " " << ec.sid << " " << m5soff(*m5) << " " << ntstr << "\n";
    			//std::cout << "--------normalize_gaps end.--------" << "\n";
    			//std::cout << "--------slide_window start.--------" << "\n";
    			//****加入的参数:slide window length
    			const int k = ctd->rco.window_length;
    			//****加入的参数:identity threshold
    			const double identity_threshold = ctd->rco.identity_threshold;
				//slide_window(nqstr, ntstr, k);
				//修改的地方####################
				slide_window2(nqstr, ntstr, newnqstr2, k, identity_threshold);
				//std::cout << "out of slide_window2" << "\n";
				//std::cout << "newnqstr:" << newnqstr << "\n";
				//std::cout << "duan point####" << "\n";
				//std::cout << "--------slide_window end.--------" << "\n";
				//meap_add_one_aln(nqstr, ntstr, m5soff(*m5), cns_table, tstr.data());//根据nqstr,ntstr的比对情况，更新cns_table中某位点的match，insert，delete, skip数量，进行了修改
				
				meap_add_one_aln(nqstr, ntstr, m5soff(*m5), cns_table, tstr.data());//修改，正常计算，不带N
				//std::cout << "meap_add_one_aln(nqstr, ntstr, m5soff(*m5), cns_table, tstr.data()); end" << "\n";
				meap_add_one_aln(newnqstr2, ntstr, m5soff(*m5), cns_table_new, tstr.data());//修改，带N
				//std::cout << "meap_add_one_aln(newnqstr2, ntstr, m5soff(*m5), cns_table_new, tstr.data()); end" << "\n";

				//修改，加入的
				cns_vec_new.add_aln(m5soff(*m5), m5send(*m5), newnqstr2, ntstr);
				//std::cout << "cns_vec_new.add_aln(m5soff(*m5), m5send(*m5), newnqstr2, ntstr); end" << "\n";

				cns_vec.add_aln(m5soff(*m5), m5send(*m5), nqstr, ntstr);//赋值a.start/end/qaln/saln
				//std::cout << "cns_vec.add_aln(m5soff(*m5), m5send(*m5), nqstr, ntstr); end" << "\n";
				
			}
		}
	}
	//****加入的参数:min_modify_coverage
    const double min_modify_coverage = ctd->rco.min_modify_coverage;
	std::vector<MappingRange> mranges, eranges;
	//std::cout << "std::vector<MappingRange> mranges, eranges; end" << "\n";
	eranges.push_back(MappingRange(0, read_size));
	//std::cout << "eranges.push_back(MappingRange(0, read_size)); end" << "\n";
	consensus_worker(cns_table, cns_table_new, ctd->id_list, cns_vec, cns_vec_new, nqstr, ntstr, eranges, ctd->rco.min_cov, ctd->rco.min_size, read_id,  cns_results, min_modify_coverage);

	//可以输出一下cns_results的元素内容，应该是一个数组

	//std::cout << "consensus_worker(cns_table, cns_table_new, ctd->id_list, cns_vec, cns_vec_new, nqstr, ntstr, eranges, ctd->rco.min_cov, ctd->rco.min_size, read_id,  cns_results); end" << "\n";
	//cns_table:存放每个位点的m,i,d,s,以及base
	//ctd->id_list
	//cns_vec:存放start/end/nqstr/ntstr
	//std::cout<<"--------consensus_one_read_can_nanopore::end--------"<<"\n";
}

} // namespace ns_meap_cns {