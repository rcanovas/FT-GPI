## FT-GPI

As presented at the TMPred web-service (http://www.ch.embnet.org/software/TMPRED_form.html): "The TMPred (Transmembrane Predictor) makes a prediction of membrane-spanning regions and their orientation. The algorithm is based on the statistical analysis of TMbase, a database of naturally occurring transmembrane proteins. The prediction is made using a combination of several weight-matrices for scoring''. The FT-GPI project offers a C++11 updated off-line version of TMPred for fasta files. Also, our version allows to modified all the internal parameters use by TMPred to detect transmembranes. For more information of how TMPred works please refer to http://www.ch.embnet.org/software/tmbase/TMBASE_doc.html


## Requirements

In order to successfully compile FT-GPI in Ubuntu-Linux machine you need: 
- An updated g++ version compiler
- An installed updated version of the cmake build system (version ≥ 3.5)
- An installed updated version of the boost library
- An installed updated version of gnuplot


##Compile

To be able to compile the codes: 
- Check that the boost library is correctly pointed in the CMakeLists.txt file
- Create the build folder if it does not exists
	- mkdir build
- Access the build folder and tun
	- cmake ..
	- make

## Methods

-[FT-GPI] :

        Usage: ./FT-GPI file_name <opt>
	  <file_name>: Name of the fasta file to be analyzed
	  
        Options:
	        -h [ --help ]                   print usage message
  		-o [ --output ] arg             pathname for output. Default: file_name.gpi
  		-t [ --training ] arg           pathname of the training set of GPI proteins.
                                  		Default: none
  		-u [ --as-table ]               if set, prints results as a table rather than
                                  		as a list. Default: not set
  		-p [ --as-percentages ]         if set, the start and end parameters are 
                                  		considered as percentages over the size of 
        		                        the sequence. Default: not set -> fixed 
                                  		integer distances
  		-b [ --both ]                   if set, checks proteins alternative model if 
                                  		strong model fails. Default: program only 
                                  		checks strong models
  		-d [ --start ] arg              maximum distance between the Nterm with the 
                                  		start of the sequence. Default: 15 amino 
                                  		acids (or 5% if p is set)
  		-e [ --end ] arg                maximum distance between the Cterm with the 
                                  		end of the sequence. Default: 10 amino acids 
                                  		(or 1% if p is set)
  		-n [ --n-score ] arg            minimum score of the Nterm. Default: 1000
  		-c [ --c-score ] arg            minimum score of the Cterm. Default: 1800
  		-i [ --i-score ] arg            maximum score allowed for all internal 
                                  		helixes. Default: -1 (not seted)
  		-k [ --i-score-percentage ] arg maximum (excluded) allowed score for internal
                                  		helixes computed as k*max-score(Nterm, 
                                  		Cterm). Default: 1.0

		TMPred Options:
  		-m [ --min-len ] arg            minimal length of transmembrane sequence. 
                                  		Default: 17
  		-M [ --max-len ] arg            maximal length of transmembrane sequence. 
                                  		Default: 35
  		-l [ --low-osl ] arg            low orientational significance level. 
                                  		Default: 80
  		-g [ --high-osl ] arg           high orientational significance level. 
                                  		Default: 200
  		-s [ --tm-osl ] arg             TM-existence significance level. Default: 500
  		-a [ --avg-osl ] arg            average orientation significance level. 
                                  		Default: 80

        
This project includes an example file in the example folder. The -g option generates the scores graphs using the gnu-iostream.h code which 
is the C++ interface to gnuplot implemented by Daniel Stahlke (http://www.stahlke.org/dan/gnuplot-iostream).

			
Note: While our tool was implemented in a Linux machine, it had been successfully tested using Cygwin (https://www.cygwin.com/) in a Window machine 
(there could be problems with the gnuplot-iostream.h code, which we recommend to comment from the tmpred.h code, including any call to its methods).

