/*
* Disktest
* Copyright (c) International Business Machines Corp., 2001
*
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*
*  Please send e-mail to yardleyb@us.ibm.com if you have
*  questions or comments.
*
*  Project Website:  TBD
*
*
* getopt.h,v 1.1 2007/07/17 12:31:52
* $Log: getopt.h,v $
* Revision 1.1  2007/07/17 12:31:52
* *** empty log message ***
*
* Revision 1.4  2005/05/04 17:54:00  mridge
* Update to version 1.2.8
*
* Revision 1.3  2002/03/30 01:32:14  yardleyb
* Major Changes:
*
* Added Dumping routines for
* data miscompares,
*
* Updated performance output
* based on command line.  Gave
* one decimal in MB/s output.
*
* Rewrote -pL IO routine to show
* correct stats.  Now show pass count
* when using -C.
*
* Minor Changes:
*
* Code cleanup to remove the plethera
* if #ifdef for windows/unix functional
* differences.
*
* Revision 1.2  2002/02/21 21:32:19  yardleyb
* Added more unix compatability
* ifdef'd function out when
* compiling for unix env. that
* have getopt
*
* Revision 1.1  2001/12/04 18:57:36  yardleyb
* This source add for windows compatability only.
*
*/
//#ifdef _WIN32
  int getopt(int argc, char** argv, const char* pszValidOpts);
//#endif /* defined WINDOWS */
