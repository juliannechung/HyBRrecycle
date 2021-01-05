  ***************************************************************************
  * All the software  contained in this library  is protected by copyright. *
  * Permission  to use, copy, modify, and  distribute this software for any *
  * purpose without fee is hereby granted, provided that this entire notice *
  * is included  in all copies  of any software which is or includes a copy *
  * or modification  of this software  and in all copies  of the supporting *
  * documentation for such software.                                        *
  ***************************************************************************
  * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED *
  * WARRANTY. IN NO EVENT, NEITHER  THE AUTHORS, NOR THE PUBLISHER, NOR ANY *
  * MEMBER  OF THE EDITORIAL BOARD OF  THE JOURNAL  "NUMERICAL ALGORITHMS", *
  * NOR ITS EDITOR-IN-CHIEF, BE  LIABLE FOR ANY ERROR  IN THE SOFTWARE, ANY *
  * MISUSE  OF IT  OR ANY DAMAGE ARISING OUT OF ITS USE. THE ENTIRE RISK OF *
  * USING THE SOFTWARE LIES WITH THE PARTY DOING SO.                        *
  ***************************************************************************
  * ANY USE  OF THE SOFTWARE  CONSTITUTES  ACCEPTANCE  OF THE TERMS  OF THE *
  * ABOVE STATEMENT AND OF THE ACCOMPANYING FILE LICENSE.txt.               *
  ***************************************************************************

   AUTHORS:

        Julianne Chung, Department of Mathematics, Virginia Tech

        Eric de Sturler, Department of Mathematics, Virginia Tech
       
        Jiahua Jiang, Department of Computer Science, Shanghai Tech
   
   REFERENCE:

       "Hybrid Projection Methods with Recycling for Inverse Problems". 
            SISC 2020.

   SOFTWARE LANGUAGE:

       MATLAB 9.9 (R2020ba)


=====================================================================
SOFTWARE
=====================================================================
The DEMO codes require the following packages:

   (1) IRTools by James Nagy, Sivia Gazzola, and Per Christian Hansen
             https://github.com/jnagy1/IRtools

   (2) AIR Tools II package from: https://github.com/jakobsj/AIRToolsII


First run startup_Recycle.m for setting paths

DEMOs on how to use HyBRrecycle

 DEMO_Deblurring.m     Sets up and runs a 2D image deblurring problem
                           corresponding to the results in Section 4.1 
                           of the paper

 DEMO_Tomo.m           Sets up and runs a 2D streaming tomgraphy problem
                           with two sets of data.  
                           This example with n = 1024 corresponds to 
                           Case 1 of Section 4.2.1 in the paper.

Supporting MATLAB functions

 HyBRrecycle.m           Main code to run HyBR with recycling
                         Syntax: [x_out, output, trunc_mats] =
                         HyBRrecycle(A, b, P, options, trunc_options, trunc_mats)
 
 recyclingGKB.m          Code to perform one step of the recycling GKB

 compression.m           Code to perform basis compression. Here we
                         provide four compression approaches: 
                         SVD, solution, sparse and RBD. 



To obtain full functionality of these codes, it is recommended to also
install 
  (3) SpaRSA by Stephen Wright, Robert Nowak, and Mario Figueiredo
             https://www.lx.it.pt/~mtf/SpaRSA/
  (4) RBD by Yanlai Chen 
             http://yanlaichen.reawritingmath.com
