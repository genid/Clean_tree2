Clean tree 2.0

Authors: Arwin Ralf, Mannis van Oven, Diego Montiel González, Sharon Wooten, Robert Lagace, Peter de Knijff, and Manfred Kayser

Department of Genetic Identification, Erasmus MC University Medical Centre Rotterdam, Rotterdam,
The Netherlands

# Requirements

    Operating system: Linux only. Tested on Ubuntu 16.04LTS, but should also work on newer version of Ubuntu. It should be easy to made it work on other Linux distributions. 
    Python, wget
    Internet connection during installation (for downloading and extracting hg19 reference genome).
    Data storage: For installation we recommend a storage capacity of > 10 GB. 

# Installation

1. Install dependencies, you can skip this step if these packages are already installed on your system
            
	    apt-get install python3.6             
            apt-get install mawk 
            apt-get install p7zip-full 
            apt-get install bwa

	-SAMtools
            We recommend the newests versions of SAMtools (e.g. > 1.4.1)

            1. wget https://github.com/samtools/samtools/releases/download/1.4.1/	samtools-1.4.1.tar.bz2 -O samtools.tar.bz2
            2. tar -xjvf samtools.tar.bz2 
            3. cd samtools-1.4.1/
            4. ./configure
            5. make
            6. make install

### Usage example for Clean Tree (BAM format)
    
        python clean_tree.py -bam file.bam -out out -r 1 -q 20 -b 95

### Usage for haplogroup prediction

	python predict_haplogroup.py -input Output_files/ -out output.hg
	
4. See complete manual at the website:
    https://www.erasmusmc.nl/genetic_identification/resources/

5. Bug report

Please email me at d.montielgonzalez@erasmusmc.nl when there is a problem getting the software up and running.
