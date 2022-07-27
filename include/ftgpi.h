/* FT-GPI 
 * Copyright (C) 2022 Rodrigo Canovas
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

#ifndef FDGPI_GPIPRED_H
#define FDGPI_GPIPRED_H


#include "tmpred.h"


namespace fdgpi {

    template<class d_type> //int if we use fixed distances, double if we use percentage distances
    class ftgpi {


    private:
        int score_nterm, score_cterm, score_iterm;
        d_type max_d_nterm, max_d_cterm;
        double k;  //Maximum score allowed on the middle helixes as percentage of the extreme scores
        bool mode;  //0 if we use fixed distances, 1 if we use percentages
        bool both_models, as_table;
        bool atg;
        uint8_t atg_d;
        std::ofstream out, out_1;
        int gpi_count, gpi_count_2;

    public:

        //Refer to FT_GPI.cpp to get a description of the parameters
        ftgpi(std::string in_name, std::string out_name, std::string training_name,
                int min_ts, int max_ts, int t1, int t2, int t3, int t4,
                double d, double e, int n, int c, int mid, double m, bool type, bool b, bool u) {
            score_nterm = n;
            score_cterm = c;
						score_iterm = mid;
            max_d_nterm = (d_type)d;
            max_d_cterm = (d_type)e;
            k = m;
            gpi_count = gpi_count_2 = 0;
            mode = type;
            both_models = b;
            as_table = u;
            if(training_name != "") {//process train data
                std::cout << "Using a training file" << std::endl;
                process_train_data(training_name, min_ts, max_ts, t1, t2, t3, t4);
            }
            //check for gpi in the input data
            gpi_predict(in_name, out_name, min_ts, max_ts, t1, t2, t3, t4);
        }

    private:

        void process_train_data(std::string f_name,
                                int m, int M, int t1, int t2, int t3, int t4) {
            tmpred tp(f_name, f_name, m, M, t1, t2, t3, t4);
            int io_score, oi_score;
            //define naive values
            max_d_nterm = max_d_cterm = 0;
            score_nterm = score_cterm = 100000; // maximum not possible score
            k = 0.5; //lower minimal reference percentage
            while(tp.process_next_sequence()) { //get data
                io_score = model_score(tp.io_helix);
                oi_score = model_score(tp.oi_helix);
                if (io_score > oi_score)
                    process_model(tp.io_helix, tp.seqlen, tp.tseq);
                else
                    process_model(tp.oi_helix, tp.seqlen, tp.tseq);
            }
        }

        void process_model(std::vector<helixtype> &model, int &seq_len, std::string &tseq) {
            int n =  model.size();
            int s_nterm, s_cterm, score_middle, max_score;
            d_type d_n, d_c;
            double per = 0.0;
            if (n > 1) {  // a valid model should have at least 2 helixes
                s_nterm = model[0].score;
                s_cterm = model[n - 1].score;
                max_score = (int)std::max(s_cterm,s_nterm);
                for (int i = 1; i < n - 1; i++) { //check that the middle section gives a valid model
                    score_middle = model[i].score;
										if (score_middle > max_score) //not valid model (check this condition)
											return;
                    per = std::max(per, score_middle * 1.0 / max_score);
                }
                d_n = compute_init_dist(model[0].nterm.position, seq_len, tseq);
                d_c = seq_len - model[n - 1].cterm.position;
                if (mode)
                    d_c = std::ceil(100 * (d_c * 1.0 / seq_len));
                if (!check_dist_cond(d_n, d_c))
                    return;
                //now I can assign rules
                k = std::max(k, per);
                max_d_nterm = std::max(max_d_nterm, d_n);
                score_nterm = (int) std::min(score_nterm, s_nterm);
                //last helix
                max_d_cterm = std::max(max_d_cterm, d_c);
                score_cterm = (int) std::min(score_cterm, s_cterm);
            }
        }

        void gpi_predict(std::string in_name, std::string out_name,
                         int m, int M, int t1, int t2, int t3, int t4) {
            tmpred tp(in_name, out_name, m, M, t1, t2, t3, t4);
            int io_score, oi_score;
            out_name += ".gpi";
            open_w(out_name, out);
            std::cout << "Creating " << out_name << std::endl;
            print_header(out);
            out_name += "_err";
            std::cout << "Creating " << out_name << std::endl;
            open_w(out_name, out_1);
            print_header(out_1, true);
            std::cout << "Processing data" << std::endl;
            while (tp.process_next_sequence()) { //get data
                atg = false;
                atg_d = 0;
                /*checks only the stronger model*/
                io_score = model_score(tp.io_helix);
                oi_score = model_score(tp.oi_helix);
                if (io_score > oi_score) {
                    if(!check_model(tp.io_helix, tp.seqlen, tp.tseq, tp.header, "str") and both_models)
                        check_model(tp.oi_helix, tp.seqlen, tp.tseq, tp.header, "alt");

                }
                else {
                    if(!check_model(tp.oi_helix, tp.seqlen, tp.tseq, tp.header, "str") and both_models)
                        check_model(tp.io_helix, tp.seqlen, tp.tseq, tp.header, "alt");
                }
            }
            std::cout << "Finish" << std::endl;
        }

        bool check_model(std::vector<helixtype> &model, int seq_len, std::string tseq,
                         std::string header, std::string m_type) {
            int n =  model.size();
            int tscore = 0;
            int n_d = 1, n_score = 1, c_d = 1, c_score = 1, m_part = 1;
            int s_nterm, s_cterm, score_middle, max_s;
            int tests = 0;
            d_type d_n, d_c, per;
            if (n < 2)//needs to have at least two helixes
                return false;
            else {
                s_nterm = model[0].score;
                s_cterm = model[n - 1].score;
                max_s = (int)std::max(s_cterm, s_nterm);
                if (n > 2) { //check if the middle part gives a valid model
                    for (int i = 1; i < n - 1; i++) {
                        score_middle = model[i].score;
												if (score_middle > max_s) {
													m_part = 0;
													break;
												}
												else {
													if (score_iterm > -1) {
														if (score_middle > score_iterm) {
															m_part = 0;
															break;
														}
													}
													else {
														per = score_middle * 1.0 / max_s;
														if (per > k) { //not valid model (need validity check of this condition)
															m_part = 0;
															break;
														}
													}
												}
										}
                } 
                d_n = compute_init_dist(model[0].nterm.position, seq_len, tseq, max_d_nterm);
                if (d_n > max_d_nterm) n_d = 0;
                if (s_nterm < score_nterm) n_score = 0;
                d_c = seq_len - model[n - 1].cterm.position;
                if (mode)
                    d_c = std::floor(100 * (d_c * 1.0 / seq_len));
                if (d_c > max_d_cterm) c_d = 0;
                if (s_cterm < score_cterm) c_score = 0;
            }
            tests = n_d + n_score + c_d + c_score + m_part;
            if (tests ==  5) {
                gpi_count++;
                if (as_table) {
                    out << std::left << std::setw(6) << std::setfill(' ') << gpi_count;
                    out << std::left << std::setw(22) << std::setfill(' ') << header.substr(1, 19);
                    out << std::left << std::setw(8) << std::setfill(' ') << seq_len;
                    out << std::left << std::setw(5) << std::setfill(' ') << (int)atg_d;
                    //style
                    if (model[0].nt_in)
                        out << std::left << std::setw(7) << std::setfill(' ') << "io";
                    else
                        out << std::left << std::setw(7) << std::setfill(' ') << "oi";
                    out << std::left << std::setw(6) << std::setfill(' ') << m_type;
                    print_model(out, model, n, tscore, as_table);
                }
                else {
                    out << gpi_count << "." << std::endl;
                    if (model[0].nt_in)
                        out << header.substr(1, 19) << " (" << m_type << ")" << std::endl << "  io model  ";
                    else
                        out << header.substr(1, 19) << " (" << m_type << ")" << std::endl << "  oi model  ";
                    out << "sequence length: " << seq_len;
                    if (atg)
                        out << "  " << (int) atg_d << " atg shift";
                    out << std::endl;
                    print_model(out, model, n, tscore);
                }
                return true;
            }
            else if (tests == 4) { //only one error allowed
                bool alt_pred = false;
                int error = 0;
                int mmid = 0; //max middle score
                if (n_d == 0) {  //Nterm distance
                    if (d_n < 1.5 * max_d_nterm) { //can change this condition
                        alt_pred = true;
                        error = 1;
                    }
                }
                else if (c_d == 0) { //Cterm distance
                    if (d_c < 1.5 * max_d_cterm) { //can change this condition
                        alt_pred = true;
                        error = 2;
                    }
                }
                else if (n_score == 0) { //Nterm score
                    if (s_nterm > 0.9 * score_nterm)  {
                        alt_pred = true;
                        error = 3;
                    }
                }
                else if (c_score == 0) { //Cterm score
                    if (s_cterm > 0.9 * score_cterm)  {
                        alt_pred = true;
                        error = 4;
                    }
                }
                else {  //Middle TM score(s)
                    for(int i=1; i< n-1; i++)
                        mmid =  std::max(mmid, model[i].score);
										if (score_iterm > -1) {
											if (mmid <= max_s and mmid <= score_iterm * 1.1) {  //can change this condition
												 alt_pred = true;
												 error = 5;
											}
										}
										else {
											per = mmid * 1.0 / max_s;
											if (mmid <= max_s and per <= (k + 0.05) ) {   //can change this condition
												alt_pred = true;
												error = 5;
											}
										}
                }
                if (alt_pred) {
                    gpi_count_2++;
                    if (as_table) {
                        out_1 << std::left << std::setw(6) << std::setfill(' ') << gpi_count_2;
                        out_1 << std::left << std::setw(8) << std::setfill(' ') << error;
                        out_1 << std::left << std::setw(22) << std::setfill(' ') << header.substr(1, 19);
                        out_1 << std::left << std::setw(8) << std::setfill(' ') << seq_len;
                        out_1 << std::left << std::setw(5) << std::setfill(' ') << (int) atg_d;
                        if (model[0].nt_in)
                            out_1 << std::left << std::setw(7) << std::setfill(' ') << "io";
                        else
                            out_1 << std::left << std::setw(7) << std::setfill(' ') << "oi";
                        out_1 << std::left << std::setw(6) << std::setfill(' ') << m_type;
                        print_model(out_1, model, n, tscore, as_table);
                    }
                    else {
                        out_1 << gpi_count_2 << "." << std::endl;
                        if (model[0].nt_in)
                            out_1 << header.substr(1, 19) << " (" << m_type << ")" << std::endl << "  io model  ";
                        else
                            out_1 << header.substr(1, 19) << " (" << m_type << ")" << std::endl << "  oi model  ";
                        out_1 << "sequence length: " << seq_len;
                        if (atg)
                            out_1 << "  " << (int) atg_d << " atg shift";
                        out_1 << "  error: " << error << std::endl;
                        print_model(out_1, model, n, tscore);
                    }
                }
            }
            return false;
        }

        void
        print_header(std::ofstream &f_out, bool err = false) {
            f_out << "GPI predicted using: (max distance, min score, term)" << std::endl;
            f_out << "(" << max_d_nterm;
            if (mode)
                f_out << "%";
            f_out << ", " << score_nterm << ", Nterm)" << std::endl;
            f_out << "(" << max_d_cterm;
            if (mode)
                f_out << "%";
            f_out << ", " << score_cterm << ", Cterm)" << std::endl;
						if (score_iterm > -1)
							f_out << score_iterm << " maximum score allowed for Middle Helixes" << std::endl << std::endl;
						else
							f_out << k << "% of the max(n-term, c-term) for Middle Helixes" << std::endl << std::endl;

            if (err) {
                f_out << "Error type: " << std::endl;
                f_out << "          1 - Start of Nterm (lower than 1.5 * limit)" << std::endl;
                f_out << "          2 - End of Cterm (lower than 1.5 * limit)" << std::endl;
                f_out << "          3 - Score of Nterm (higher than 0.9 * limit)" << std::endl;
                f_out << "          4 - Score of Cterm (higher than 0.9 * limit)" << std::endl;
                f_out << "          5 - Score Middle TM(s) (add 5% to the k limit or add 10% to the score limit)" << std::endl << std::endl;
            }

            if (as_table) {
                f_out << std::left << std::setw(6) << std::setfill(' ') << " ";
                if (err)
                    f_out << std::left << std::setw(8) << std::setfill(' ') << "Error";
                f_out << std::left << std::setw(22) << std::setfill(' ') << "Protein name";
                f_out << std::left << std::setw(8) << std::setfill(' ') << "length";
                f_out << std::left << std::setw(5) << std::setfill(' ') << "atg";
                f_out << std::left << std::setw(7) << std::setfill(' ') << "model";
                f_out << std::left << std::setw(6) << std::setfill(' ') << "type";
                f_out << std::left << std::setw(10) << std::setfill(' ') << "t. score";
                f_out << std::left << std::setw(21) << std::setfill(' ') << "Nterm";
                f_out << std::left << std::setw(21) << std::setfill(' ') << "Cterm";
                f_out << std::left << std::setw(21) << std::setfill(' ') << "Middle -->";
                f_out << std::endl;
                if (err)
                    f_out << std::left << std::setw(68) << std::setfill(' ') << " ";
                else
                    f_out << std::left << std::setw(60) << std::setfill(' ') << " ";
                f_out << std::right << std::setw(21) << std::setfill(' ') << "(From, To, Score)";
                f_out << std::right << std::setw(21) << std::setfill(' ') << "(From, To, Score)";
                f_out << std::right << std::setw(21) << std::setfill(' ') << "(From, To, Score)";
                f_out << std::endl;
            }
        }

        d_type
        compute_init_dist(int n_pos, int &seq_len, std::string &seq, d_type limit = 0) {
            d_type d = n_pos;
            if (limit == 0) {
                if (mode)   limit = 10;
                else        limit = 30;
            }
            int spos = 0;
            if (mode)  //we are using percentage
                d = std::floor(100 * d * 1.0 / seq_len);
            while (d > limit) {
                atg = true;
                atg_d ++;
                spos++;
                while (spos <= n_pos and seq[spos] != 'M') //find next ATG start
                    spos++;
                if (spos > n_pos or spos > 30)
                    break;
                d = n_pos - spos;
                if (mode)
                    d = std::floor(100 * d * 1.0 / seq_len);
            }
            return d;
        }

        bool
        check_dist_cond(d_type d_n, d_type d_c) {
            if (mode) {
                if (d_n > 10 or d_c > 10) //we assumed that 10% is too big
                    return false;
            }
            else {
                if (d_n > 30 or d_c > 30) //we assumed that 30 is too big
                    return false;
            }
            return true;
        }

    };
}

#endif //FDGPI_GPIPRED_H
