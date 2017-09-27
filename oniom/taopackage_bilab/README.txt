==============================================================================
===                                  TAO                                   ===
==============================================================================


The Toolkit to Assist ONIOM (TAO) calculations is designed to assist in the 
different stages of an ONIOM QM/MM study of biomolecules, including input 
file preparation and checking, job monitoring, production calculations, 
and results analysis. 

TAO is an open-source toolkit in PERL.

Please check out http://www.chem.wayne.edu/schlegel/Software.html

To install this package on your UNIX/LINUX platform, please follow
this three step installation:

1)
Move this toolkit package to a location you usually install  
application softwares.


2)
Edit install.sh in the home folder of TAO (taopackage).
Change path ~/bin in line 14

USERPATH=~/bin

to the path that you want the symbolic links to TAO
scripts to be installed. e.g.

USERPATH=/home/myhome/bin/oniomtool

Please make sure this path is in your search PATH.

3)
Run ./install.sh in the home folder of TAO (taopackage) to 
install this package.


Notice: If you use tcsh, you may need to log out and log back into your account 
        right after you run ./install.sh. The toolkit will run afterward.


Please contact Dr. Peng Tao (tao.21@osu.edu) for any further questions.

The author give special thanks Brian T. Psciuk.


Copyright (c) 2007-2009 Peng Tao


