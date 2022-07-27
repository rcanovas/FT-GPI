/* TMPred - Transmembrane Region Predictor
 * Copyright (C) 2018 Rodrigo Canovas
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses/ .
 * */


#ifndef FDGPI_TMPRED_H
#define FDGPI_TMPRED_H


#include "tables.h"
#include "ureadseq.h"
#include <iomanip>

namespace fdgpi {

    struct featuretype {
        int  position;
        int score;
    };

    struct helixtype {
        featuretype center, nterm, cterm, sh_cterm, sh_nterm;
        int score;
        bool nt_in;
    };

    int model_score(std::vector<helixtype> &model) {
        int sum = 0, count = model.size();
        for (int i = 0; i < count; i++)
            sum += model[i].score;
        return sum;
    }

    void print_model(std::ofstream &f, std::vector<helixtype> &model, int &count, int &totscore, bool as_table = false) {
        int n = model.size();
        if (totscore == 0)
            totscore = model_score(model);
        if (as_table) {
            f << std::left << std::setw(10) << std::setfill(' ') << totscore; //t. score
            //NTerm
            f << std::left << std::setw(7) << model[0].nterm.position; //from
            f << std::left << std::setw(7) << model[0].cterm.position; //to
            f << std::left << std::setw(7) << model[0].score; //score
            //CTerm
            f << std::left << std::setw(7) << model[n - 1].nterm.position; //from
            f << std::left << std::setw(7) << model[n - 1].cterm.position; //to
            f << std::left << std::setw(7) << model[n - 1].score; //score
            //Middle
            for (int i = 1; i < (n - 1); i++) {
                f << std::left << std::setw(7) << model[i].nterm.position; //from
                f << std::left << std::setw(7) << model[i].cterm.position; //to
                f << std::left << std::setw(7) << model[i].score; //score
            }
        }
        else {
            f << "  " << count << " strong transmembrane helices, total score : "
              << totscore << std::endl;
            f << "  #  from   to  length   score  orientation " << std::endl;
            for (int i = 1; i <= count; i++) {
                f << "  " << i << "   " << model[i - 1].nterm.position << "   "
                  << model[i - 1].cterm.position;
                f << "  (" << model[i - 1].cterm.position - model[i - 1].nterm.position + 1
                  << ")      " << model[i - 1].score;
                if (model[i - 1].nt_in)
                    f << "    i-o" << std::endl;
                else
                    f << "    o-i" << std::endl;
            }
        }
        f << std::endl;
    }


    class tmpred {

    public:
        std::vector<helixtype> io_helix, oi_helix; //used to store the helixes and the models at the end
        int seqlen;
        std::string tseq;
        std::string header;

    private:

        int thres1;      //low orientational significance level
        int thres2;      //high orientational significance level
        int thres3;      //TM-existence significance level
        int thres4;      //average orientation significance level


        std::vector<int > iom, ion, ioc, oim, oin, oic;
        std::vector<int > io_score, oi_score;

        std::vector<double> scale; //frequency vector used to scale the tables

        int minwidth, maxwidth; //minimal and maximal length of transmembrane sequence
        int io_count, oi_count, maxhelix = 75;
        int seq_count;

        std::ofstream out;
        std::ifstream f_in;
        rwSeq rws; //used for reading sequences
        bool first;

    public:
        tmpred(std::string in_name, std::string out_name,
               int min_ts, int max_ts, int t1, int t2, int t3, int t4) {
            first = true;
            thres1 = t1;
            thres2 = t2;
            thres3 = t3;
            thres4 = t4;
            minwidth = min_ts;
            maxwidth = max_ts;
            {
                std::ifstream in_tmp(in_name, std::ios::in | std::ios::binary);
                if (!in_tmp) {
                    std::cerr << "Failed to open file " << in_name << std::endl;
                    exit(1);
                }
                f_in.swap(in_tmp);
            }
        }

        void process_all() {
            while(process_next_sequence()){ }
        }

        //if tm_tp ==false then the transmembrane_topology step is skipped
        bool process_next_sequence() {
            seqlen = MAXSEQLEN;
            if (rws.readNextSeq(f_in, seqlen, tseq, header, first)) {
                restart_variable();
                seq_count++;
                for (auto & c: tseq) c = std::toupper(c);
                compute_tmpred();
                return true;
            }
            else
                return false;
        }

    private:

        //Computes all the tmpred information of the current tseq
        void compute_tmpred() {
            compute_tables();  //creates the profile arrays from the tables
            compute_scores();  //computes io/oi_scores
            predict();   //predict helixes
            transmembrane_topology();
        }

        void restart_variable() {
            iom.clear();
            ion.clear();
            ioc.clear();
            oim.clear();
            oin.clear();
            oic.clear();
            io_score.clear();
            oi_score.clear();
            io_helix.clear();
            oi_helix.clear();
        }

        void compute_tables() {
            std::vector<double> scale; //frequency vector used to scale the tables
            scale = aa_freq;  //can be change for aa_freq_old or aa_freq_tm
            std::vector<std::vector<double> > m_io_ce, m_oi_ce, m_io_nt, m_oi_nt, m_io_ct, m_oi_ct;
            scaletables(scale, m_io_ce, m_oi_ce, m_io_nt, m_oi_nt, m_io_ct, m_oi_ct);
            make_profile(m_io_ce, 11, iom);
            make_profile(m_io_nt, 6, ion);
            make_profile(m_io_ct, 16, ioc);
            make_profile(m_oi_ce, 11, oim);
            make_profile(m_oi_nt, 6, oin);
            make_profile(m_oi_ct, 16, oic);
        }

        void compute_scores() {
            make_curve(true);
            make_curve(false);
        }

        void predict() {
            int start = 1;
            bool found = false;
            helixtype helix;
            io_count = 0;

            while (start < seqlen and io_count <= (maxhelix / 2)) {
                found = find_helix(start, helix, true);
                if (found) {
                    helix.nt_in = true;
                    io_count ++;
                    io_helix.push_back(helix);
                }
            }
            start = 1;
            oi_count = 0;
            while (start < seqlen and oi_count <= (maxhelix / 2)) {
                found = find_helix(start, helix, false);
                if (found) {
                    helix.nt_in = false;
                    oi_count++;
                    oi_helix.push_back(helix);
                }
            }
        }


        //Returns the maximum value in the array prf and its position in the parameter p
        int64_t  findmax(std::vector<int> &prf, int from, int to, int &p) {
            int m = prf[from]; p = from;
            if (to < from)
                return m;
            if (to > prf.size())
                to =  prf.size() - 1;
            for (int i = from + 1; i <= to; i++) {
                if (prf[i] > m) {
                    m = prf[i];
                    p = i;
                }
            }
            return m;
        }

        // Creates the profile of a table and write it in p
        // the refpos is the position to be labeled
        void make_profile(std::vector<std::vector<double>> d, int refpos, std::vector<int> &p) {
            std::vector<int> profile;
            profile.push_back(0);
            uint8_t ncol = n_col; //21
            int pos, rc;
            double m;
            for (int k = 1; k <= seqlen; k ++) {
                m = 0;
                for (int i = 1; i <= ncol; i ++) {
                    pos = k + i - refpos; //+2
                    if ( (pos >= 1) and (pos <= seqlen)) {
                        rc = rowcode.find(tseq[pos - 1]);
                        if (rc == std::string::npos)
                            m += d[0][i];
                        else
                            m += d[rowcode.find(tseq[pos - 1]) + 1][i];
                    }
                }
                profile.push_back(round(m * 100));
            }
            profile.swap(p);
        }

        void make_curve(bool io) {
            std::vector<int> *s, *m, *n, *c;
            int value;
            if (io) {
                s = &io_score; m = &iom; n = &ion; c = &ioc;
            }
            else {
                s = &oi_score; m = &oim; n = &oin; c = &oic;
            }
            int dummy;
            int minhw = minwidth / 2, maxhw = (maxwidth + 1) / 2;
            s->push_back(0); // first value equal to zero
            for (int i = 1; i <= seqlen; i++) {
                if (i <= minhw  or i > (seqlen - minhw))
                    s->push_back(0);
                else {
                    value = (*m)[i];
                    value += findmax((*n), std::max(i - maxhw, 1), i - minhw, dummy);
                    value += findmax((*c), i + minhw, std::min(i + maxhw, seqlen), dummy);
                    s->push_back(value);
                }
            }
        }

        //Search for an helix in tseq starting from start.
        // io indicates if we are using io or oi parameters
        bool find_helix(int &start, helixtype &helix, bool io) {
            int minhw = minwidth / 2, maxhw = (maxwidth + 1) / 2;
            bool found = false, done;
            int i = start, max_value, dummy, j;
            std::vector<int> *s, *m, *n, *c;
            if (io) {
                s = &io_score; m = &iom; n = &ion; c = &ioc;
            }
            else {
                s = &oi_score; m = &oim; n = &oin; c = &oic;
            }
            if (i <= minhw)
                i = minhw + 1;
            while ((i <= seqlen - minhw)  and !found) {
                max_value = findmax((*s), i - minhw, i + maxhw, dummy);
                if (((*s)[i] == max_value) and ((*s)[i] > 0)) {
                    found = true;
                    helix.center.position = i;
                    helix.center.score = (*m)[i]; //determine center of TM-segment
                    helix.nterm.score = findmax((*n),
                                                std::max(i - maxhw, 1),
                                                i - minhw, helix.nterm.position);  //optimal..
                    helix.cterm.score = findmax((*c),
                                               i + minhw,
                                               std::min(i + maxhw, seqlen), helix.cterm.position); //..termini
                    j = i - minhw; //determine nearest N-terminus
                    done = false;
                    while((j - 1 >= 1) and ((j - 1) >= i - maxhw) and !done) {
                        if ((*n)[j - 1] > (*n)[j])
                            j--;
                        else
                            done = true;
                    }
                    helix.sh_nterm.score = (*n)[j];
                    helix.sh_nterm.position = j;
                    j = i + minhw;  //determine nearest C-terminus
                    done = false;
                    while ( (j + 1 <= seqlen) and (j + 1 <= i + maxhw) and !done) {
                        if (((*c)[j + 1] > (*c)[j]))
                            j++;
                        else
                            done = true;
                    }
                    helix.sh_cterm.score = (*c)[j];
                    helix.sh_cterm.position = j;
                }
                i ++;
            }
            if (found) {
                start = helix.sh_cterm.position + 1;
                helix.score = helix.center.score + helix.nterm.score + helix.cterm.score;
            }
            else
                start = seqlen;
            return found;
        }

        //mode: the array to search in
        //ct: where to start searching, ct =0 if not found
        //maxct: no of TM-segs in array
        //ntpos: earliest pos. of N-term
        //thres: threshold for ignoring
        void findnext_helix(std::vector<helixtype> &model, int &ct, int &maxct, int &ntpos, int &thres){
            if (ct > maxct)
                ct = 0;
            else {
                while ((ct < maxct) and ((model[ct - 1].score < thres) or (model[ct - 1].sh_nterm.position < ntpos)))
                    ct++;
                if ((model[ct - 1].score < thres) or (model[ct - 1].sh_nterm.position < ntpos))
                    ct = 0;
            }
        }

        //ori_io: is N-term inside?,
        //model: output helixes
        //model_ct: number of segments in model
        void compute_model(bool ori_io, std::vector<helixtype> &model, int &model_ct) {
            int io_ct = 1, oi_ct = 1, nt_pos = 1;
            bool done = false;
            model_ct = 0;
            while (!done) {
                if (ori_io) {
                    findnext_helix(io_helix, io_ct, io_count, nt_pos, thres3);
                    if (io_ct == 0)
                        done = true;
                    else {
                        model_ct ++;
                        model.push_back(io_helix[io_ct - 1]);
                        nt_pos = model[model_ct - 1].sh_cterm.position + 1;
                    }
                }
                else {
                    findnext_helix(oi_helix, oi_ct, oi_count, nt_pos, thres3);
                    if (oi_ct == 0)
                        done = true;
                    else {
                        model_ct ++;
                        model.push_back(oi_helix[oi_ct - 1]);
                        nt_pos = model[model_ct - 1].sh_cterm.position + 1;
                    }
                }
                ori_io = !ori_io;
            }
        }

        void transmembrane_topology() {
            std::vector<helixtype> in_model, out_model;
            int in_mod_score, out_mod_score, in_mod_count, out_mod_count;
            compute_model(true, in_model, in_mod_count);
            in_mod_score = model_score(in_model);
            compute_model(false, out_model, out_mod_count);
            out_mod_score = model_score(out_model);
            io_helix.swap(in_model);
            oi_helix.swap(out_model);
        }

    };
}


#endif //FDGPI_TMPRED_H
