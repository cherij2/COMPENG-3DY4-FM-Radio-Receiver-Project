## Real-time SDR for mono/stereo FM and RDS

The main objective of the project is to navigate a complex speciﬁcation and develop an understanding of the challenges that must be addressed for a real-time implementation of a computing system operating in a form factor-constrained environment. 

The project description is available in the project [document](doc/3dy4-project-2024.pdf). The unique constraints for each group for different modes of operation (i.e., custom sample rates) are available [here](doc/3dy4-constraints-group-4.pdf).

All the project source code must be submitted before 7 p.m. on March 28. The project cross-examinations and oral presentations will run in the week of April 1, and the detailed final project report is due before 7 p.m. on April 5. 


If you want to test RDS without multithreading uncomment lines 78 to 301 in the project.cpp (Bit recovery from RRC coded, Manchester encoding not fully coded). 
Multithreading works on all modes, working completely with samples. One underrun occurs around block 35 on rare occasion. Tested on March 28, 2024 live with the PI (NO UNDERRUNS). Tried changing queue size and block size.
