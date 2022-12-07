.. :changelog:

History
-------

0.3.1 (12/7/2022)
~~~~~~~~~~~~~~~~~~

* Hotfix for multifile inputs. Output file names are now auto-generated for more than one input.

0.3.0 (12/7/2022)
~~~~~~~~~~~~~~~~~~

* Fixed reverse conversions
* Periods can now be in contig names
* New file formats:
 - .junc files from leafcutter/regtools
 - .tab files from STAR
* Now can convert multiple files at once using one vci
* Haploid contigs in diploid vcis are handled appropriately (note: these frankenstein vcis still need to be made manually)
* When reverse-converting a sam/bam from diploid coordinates back to reference haploid coordinates, reads have a 'vA' tag set to record the allele to which it originally aligned. This tag mirrors the 'vA' tag set by STAR.

0.2.9 (10/01/2019)
~~~~~~~~~~~~~~~~~~

* Fixed error of not setting correct header information in B/SAM.

0.2.8 (09/23/2019)
~~~~~~~~~~~~~~~~~~

* Fixed error of not parsing VCI file on bam/sam conversion.

0.2.5 (06/06/2018)
~~~~~~~~~~~~~~~~~~

* Fixed error of always generating index file

0.2.4 (06/05/2018)
~~~~~~~~~~~~~~~~~~

* Automtically generates file index if not found

0.2.3 (06/05/2018)
~~~~~~~~~~~~~~~~~~

* Fixed end of contig tree mapping

0.2.2 (06/04/2018)
~~~~~~~~~~~~~~~~~~

* Fixed extract transcripts

0.2.1 (05/30/2018)
~~~~~~~~~~~~~~~~~~

* Hotfix to solve gtf parsing problem in python 3

0.2.0 (05/09/2018)
~~~~~~~~~~~~~~~~~~

* Support for diploid VCF
* Reimplemented to support G2G VCI (Variant Call Information)
* Support python3

0.1.31 (02/17/2017)
~~~~~~~~~~~~~~~~~~~

* Final release for custom haploid reconstruction

0.1.29 (05/17/2016)
~~~~~~~~~~~~~~~~~~~

* Fixed parsing issue with GATK-generated VCF files

0.1.27 (05/04/2016)
~~~~~~~~~~~~~~~~~~~

* Uploaded the package to Anaconda Cloud
* Fixed Travis CI fail

0.1.23 (02/15/2016)
~~~~~~~~~~~~~~~~~~~

* Setup Travis CI for automated building

0.1.22 (02/15/2016)
~~~~~~~~~~~~~~~~~~~

* Updated documentation

0.1.20 (02/03/2016)
~~~~~~~~~~~~~~~~~~~

* Released to public with documentation at g2gtools.readthedocs.org

0.1.0 (02/09/2015)
~~~~~~~~~~~~~~~~~~

* Started github repository
