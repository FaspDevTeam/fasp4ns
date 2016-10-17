#UMFPACKDIR = /Users/zhangcs/Packages
UMFPACKINCLUDE = -I$(UMFPACKDIR)/include
UMFPACKLIB = -L$(UMFPACKDIR)/lib -lumfpack -lamd -lcholmod -lcolamd -lcamd -lccolamd #-lSuiteSparse
