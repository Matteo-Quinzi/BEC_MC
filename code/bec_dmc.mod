  è=     k820309    ,          2021.5.0    Vub                                                                                                          
       bec_dmc.f90 BEC_DMC                      @                              
       %         @                                                   
       #A    #RI    #RJ              
                                      
                
                                                    
    p          p            p                                    
                                                    
 	   p          p            p                          %         @                                                   
       #A    #RIJ              
                                      
                
                                      
      %         @                               	                    
       #A 
   #XIJ    #RIJ              
                                 
     
                
                                      
                
                                      
      %         @                                                   
       #B0    #B1    #R_IDX    #R2              
                                      
                
                                      
                
                                      
                
                                      
      %         @                                                   
       #B0    #B1    #R              
                                      
                
                                      
                
                                                    
    p          p            p                          %         @                                                    
       #A    #B0    #B1    #N_AT    #COORDS              
                                      
                
                                      
                
                                      
                
                                                     
                                                     
      p        4 5 O p        p          4 5 O p          p            4 5 O p          p                          %         @                                                           #LHS    #RHS              
                                                     #MPI_INFO              
                                                     #MPI_INFO    %         @                                                            #LHS !   #RHS #             
                                  !                   #MPI_REQUEST "             
                                  #                   #MPI_REQUEST "   %         @                                $                           #LHS %   #RHS '             
                                  %                   #MPI_COMM &             
                                  '                   #MPI_COMM &   %         @                                (                           #LHS )   #RHS +             
                                  )                   #MPI_WIN *             
                                  +                   #MPI_WIN *   %         @                                ,                           #LHS -   #RHS /             
                                  -                   #MPI_FILE .             
                                  /                   #MPI_FILE .   %         @                                0                           #LHS 1   #RHS 3             
                                  1                   #MPI_MESSAGE 2             
                                  3                   #MPI_MESSAGE 2   %         @                                4                           #LHS 5   #RHS 7             
                                  5                   #MPI_ERRHANDLER 6             
                                  7                   #MPI_ERRHANDLER 6   %         @                                8                           #LHS 9   #RHS ;             
                                  9                   #MPI_GROUP :             
                                  ;                   #MPI_GROUP :   %         @                                <                           #LHS =   #RHS ?             
                                  =                   #MPI_DATATYPE >             
                                  ?                   #MPI_DATATYPE >   %         @                                @                           #LHS A   #RHS C             
                                  A                   #MPI_OP B             
                                  C                   #MPI_OP B   %         @                                D                           #LHS E   #RHS F             
                                  E                   #MPI_INFO              
                                  F                   #MPI_INFO    %         @                                G                           #LHS H   #RHS I             
                                  H                   #MPI_REQUEST "             
                                  I                   #MPI_REQUEST "   %         @                                J                           #LHS K   #RHS L             
                                  K                   #MPI_COMM &             
                                  L                   #MPI_COMM &   %         @                                M                           #LHS N   #RHS O             
                                  N                   #MPI_WIN *             
                                  O                   #MPI_WIN *   %         @                                P                           #LHS Q   #RHS R             
                                  Q                   #MPI_FILE .             
                                  R                   #MPI_FILE .   %         @                                S                           #LHS T   #RHS U             
                                  T                   #MPI_MESSAGE 2             
                                  U                   #MPI_MESSAGE 2   %         @                                V                           #LHS W   #RHS X             
                                  W                   #MPI_ERRHANDLER 6             
                                  X                   #MPI_ERRHANDLER 6   %         @                                Y                           #LHS Z   #RHS [             
                                  Z                   #MPI_GROUP :             
                                  [                   #MPI_GROUP :   %         @                                \                           #LHS ]   #RHS ^             
                                  ]                   #MPI_DATATYPE >             
                                  ^                   #MPI_DATATYPE >   %         @                                _                           #LHS `   #RHS a             
                                  `                   #MPI_OP B             
                                  a                   #MPI_OP B                                              b     
                                                   c     ACOS (        `                                d                                    
    #N_WALK e   #N_AT f   #N_DIM g       p        4 5 O p        p        4 5 O p        p          4 5 O p          4 5 O p          4 5 O p            4 5 O p          4 5 O p          4 5 O p                                    
                                 e                     
                                 f                     
                                 g           (        `                                h                    
                
    #A i   #B0 j   #B1 k   #N_AT l   #COORDS m     p        4 5 O p        p          4 5 O p          p            4 5 O p          p                                    
  @                              i     
                
  @                              j     
                
  @                              k     
                
                                 l                    
                                 m                    
      p        4 5  p        r l   p          4 5  p        r l     p            4 5  p        r l     p                          %         @                                n                    
       #A o   #B0 p   #B1 q   #N_AT r   #COORDS s   #NEW_COORDS t   #F_DRIV u   #F_DRIV_NEW v   #E_LOC w   #E_LOC_NEW x   #DT y             
  @                              o     
                
  @                              p     
                
  @                              q     
                
                                 r                    
                                 s                    
      p        4 5  p        r r   p          4 5  p        r r     p            4 5  p        r r     p                                   
                                 t                    
      p        4 5  p        r r   p          4 5  p        r r     p            4 5  p        r r     p                                   
                                 u                    
      p        4 5  p        r r   p          4 5  p        r r     p            4 5  p        r r     p                                   
                                 v                    
      p        4 5  p        r r   p          4 5  p        r r     p            4 5  p        r r     p                                    
                                 w     
                
                                 x     
                
                                 y     
      #         @                                   z                 	   #A {   #B0 |   #B1 }   #N_AT ~   #DT    #COORDS    #F    #E_LOC    #ACCEPTED              
  @                              {     
                
  @                              |     
                
  @                              }     
                
  @                              ~                     
  @                                   
               
D @                                                  
       p        4 5  p        r ~   p          4 5  p        r ~     p            4 5  p        r ~     p                                   
D @                                                  
       p        4 5  p        r ~   p          4 5  p        r ~     p            4 5  p        r ~     p                                    
D @                                   
                 
D                                             #         @                                                    	   #N_WALK    #N_MAX    #N_AT    #CONFIGURATIONS    #WALK_EN    #F_DRIV    #FLAG    #DT    #ER              
D                                                      
                                                      
                                                     
D                                                    
         p        4 5  p        r    p        4 5  p        r    p          4 5  p        r      4 5  p        r      p            4 5  p        r      4 5  p        r      p                                   
D                                                    
     p          4 5  p        r        4 5  p        r                               
D                                                    
         p        4 5  p        r    p        4 5  p        r    p          4 5  p        r      4 5  p        r      p            4 5  p        r      4 5  p        r      p                                   
D                                                          p          4 5  p        r        4 5  p        r                                
                                      
                
                                      
                       @                           B     '                    #MPI_VAL                                                                                    @                           >     '                    #MPI_VAL                                                                                    @                           :     '                    #MPI_VAL                                                                                    @                           6     '                    #MPI_VAL                                                                                    @                           2     '                    #MPI_VAL                                                                                    @                           .     '                    #MPI_VAL                                                                                    @                           *     '                    #MPI_VAL                                                                                    @                           &     '                    #MPI_VAL                                                                                    @                           "     '                    #MPI_VAL                                                                                    @                                '                    #MPI_VAL                                                                                fn#fn    ¼   @   J   BEC_VMC    ü   g       F+BEC_VMC    c  @   a   F%A+BEC_VMC    £     a   F%RI+BEC_VMC    7     a   F%RJ+BEC_VMC    Ë  `       FAST_F+BEC_VMC !   +  @   a   FAST_F%A+BEC_VMC #   k  @   a   FAST_F%RIJ+BEC_VMC )   «  i       FAST_FIRST_DER_F+BEC_VMC +     @   a   FAST_FIRST_DER_F%A+BEC_VMC -   T  @   a   FAST_FIRST_DER_F%XIJ+BEC_VMC -     @   a   FAST_FIRST_DER_F%RIJ+BEC_VMC *   Ô  s       FIRST_DERG_OVER_G+BEC_VMC -   G  @   a   FIRST_DERG_OVER_G%B0+BEC_VMC -     @   a   FIRST_DERG_OVER_G%B1+BEC_VMC 0   Ç  @   a   FIRST_DERG_OVER_G%R_IDX+BEC_VMC -     @   a   FIRST_DERG_OVER_G%R2+BEC_VMC    G  g       G+BEC_VMC    ®  @   a   G%B0+BEC_VMC    î  @   a   G%B1+BEC_VMC    .     a   G%R+BEC_VMC %   Â  }       LOCAL_ENERGY+BEC_VMC '   ?  @   a   LOCAL_ENERGY%A+BEC_VMC (     @   a   LOCAL_ENERGY%B0+BEC_VMC (   ¿  @   a   LOCAL_ENERGY%B1+BEC_VMC *   ÿ  @   a   LOCAL_ENERGY%N_AT+BEC_VMC ,   ?	  ø   a   LOCAL_ENERGY%COORDS+BEC_VMC %   7
  b       INFOEQ+MPI_CONSTANTS )   
  V   a   INFOEQ%LHS+MPI_CONSTANTS )   ï
  V   a   INFOEQ%RHS+MPI_CONSTANTS (   E  b       REQUESTEQ+MPI_CONSTANTS ,   §  Y   a   REQUESTEQ%LHS+MPI_CONSTANTS ,      Y   a   REQUESTEQ%RHS+MPI_CONSTANTS %   Y  b       COMMEQ+MPI_CONSTANTS )   »  V   a   COMMEQ%LHS+MPI_CONSTANTS )     V   a   COMMEQ%RHS+MPI_CONSTANTS $   g  b       WINEQ+MPI_CONSTANTS (   É  U   a   WINEQ%LHS+MPI_CONSTANTS (     U   a   WINEQ%RHS+MPI_CONSTANTS %   s  b       FILEEQ+MPI_CONSTANTS )   Õ  V   a   FILEEQ%LHS+MPI_CONSTANTS )   +  V   a   FILEEQ%RHS+MPI_CONSTANTS (     b       MESSAGEEQ+MPI_CONSTANTS ,   ã  Y   a   MESSAGEEQ%LHS+MPI_CONSTANTS ,   <  Y   a   MESSAGEEQ%RHS+MPI_CONSTANTS +     b       ERRHANDLEREQ+MPI_CONSTANTS /   ÷  \   a   ERRHANDLEREQ%LHS+MPI_CONSTANTS /   S  \   a   ERRHANDLEREQ%RHS+MPI_CONSTANTS &   ¯  b       GROUPEQ+MPI_CONSTANTS *     W   a   GROUPEQ%LHS+MPI_CONSTANTS *   h  W   a   GROUPEQ%RHS+MPI_CONSTANTS )   ¿  b       DATATYPEEQ+MPI_CONSTANTS -   !  Z   a   DATATYPEEQ%LHS+MPI_CONSTANTS -   {  Z   a   DATATYPEEQ%RHS+MPI_CONSTANTS #   Õ  b       OPEQ+MPI_CONSTANTS '   7  T   a   OPEQ%LHS+MPI_CONSTANTS '     T   a   OPEQ%RHS+MPI_CONSTANTS &   ß  b       INFONEQ+MPI_CONSTANTS *   A  V   a   INFONEQ%LHS+MPI_CONSTANTS *     V   a   INFONEQ%RHS+MPI_CONSTANTS )   í  b       REQUESTNEQ+MPI_CONSTANTS -   O  Y   a   REQUESTNEQ%LHS+MPI_CONSTANTS -   ¨  Y   a   REQUESTNEQ%RHS+MPI_CONSTANTS &     b       COMMNEQ+MPI_CONSTANTS *   c  V   a   COMMNEQ%LHS+MPI_CONSTANTS *   ¹  V   a   COMMNEQ%RHS+MPI_CONSTANTS %     b       WINNEQ+MPI_CONSTANTS )   q  U   a   WINNEQ%LHS+MPI_CONSTANTS )   Æ  U   a   WINNEQ%RHS+MPI_CONSTANTS &     b       FILENEQ+MPI_CONSTANTS *   }  V   a   FILENEQ%LHS+MPI_CONSTANTS *   Ó  V   a   FILENEQ%RHS+MPI_CONSTANTS )   )  b       MESSAGENEQ+MPI_CONSTANTS -     Y   a   MESSAGENEQ%LHS+MPI_CONSTANTS -   ä  Y   a   MESSAGENEQ%RHS+MPI_CONSTANTS ,   =  b       ERRHANDLERNEQ+MPI_CONSTANTS 0     \   a   ERRHANDLERNEQ%LHS+MPI_CONSTANTS 0   û  \   a   ERRHANDLERNEQ%RHS+MPI_CONSTANTS '   W  b       GROUPNEQ+MPI_CONSTANTS +   ¹  W   a   GROUPNEQ%LHS+MPI_CONSTANTS +     W   a   GROUPNEQ%RHS+MPI_CONSTANTS *   g  b       DATATYPENEQ+MPI_CONSTANTS .   É  Z   a   DATATYPENEQ%LHS+MPI_CONSTANTS .   #  Z   a   DATATYPENEQ%RHS+MPI_CONSTANTS $   }  b       OPNEQ+MPI_CONSTANTS (   ß  T   a   OPNEQ%LHS+MPI_CONSTANTS (   3  T   a   OPNEQ%RHS+MPI_CONSTANTS      @       PI    Ç  =       ACOS       ¥      GAUSSIAN_RNG $   ©!  @   a   GAUSSIAN_RNG%N_WALK "   é!  @   a   GAUSSIAN_RNG%N_AT #   )"  @   a   GAUSSIAN_RNG%N_DIM    i"  5      DRIVING_FORCE     #  @   a   DRIVING_FORCE%A !   Þ#  @   a   DRIVING_FORCE%B0 !   $  @   a   DRIVING_FORCE%B1 #   ^$  @   a   DRIVING_FORCE%N_AT %   $    a   DRIVING_FORCE%COORDS    ®%  Ë       DMC_ACC_PROB    y&  @   a   DMC_ACC_PROB%A     ¹&  @   a   DMC_ACC_PROB%B0     ù&  @   a   DMC_ACC_PROB%B1 "   9'  @   a   DMC_ACC_PROB%N_AT $   y'    a   DMC_ACC_PROB%COORDS (   (    a   DMC_ACC_PROB%NEW_COORDS $   )    a   DMC_ACC_PROB%F_DRIV (   ©*    a   DMC_ACC_PROB%F_DRIV_NEW #   ¹+  @   a   DMC_ACC_PROB%E_LOC '   ù+  @   a   DMC_ACC_PROB%E_LOC_NEW     9,  @   a   DMC_ACC_PROB%DT %   y,         ONE_WALKER_DIFFUSION '   -  @   a   ONE_WALKER_DIFFUSION%A (   V-  @   a   ONE_WALKER_DIFFUSION%B0 (   -  @   a   ONE_WALKER_DIFFUSION%B1 *   Ö-  @   a   ONE_WALKER_DIFFUSION%N_AT (   .  @   a   ONE_WALKER_DIFFUSION%DT ,   V.    a   ONE_WALKER_DIFFUSION%COORDS '   f/    a   ONE_WALKER_DIFFUSION%F +   v0  @   a   ONE_WALKER_DIFFUSION%E_LOC .   ¶0  @   a   ONE_WALKER_DIFFUSION%ACCEPTED    ö0  °       BRANCHING !   ¦1  @   a   BRANCHING%N_WALK     æ1  @   a   BRANCHING%N_MAX    &2  @   a   BRANCHING%N_AT )   f2    a   BRANCHING%CONFIGURATIONS "   ò3  ¼   a   BRANCHING%WALK_EN !   ®4    a   BRANCHING%F_DRIV    :6  ¼   a   BRANCHING%FLAG    ö6  @   a   BRANCHING%DT    67  @   a   BRANCHING%ER %   v7  ]       MPI_OP+MPI_CONSTANTS -   Ó7  H   a   MPI_OP%MPI_VAL+MPI_CONSTANTS +   8  ]       MPI_DATATYPE+MPI_CONSTANTS 3   x8  H   a   MPI_DATATYPE%MPI_VAL+MPI_CONSTANTS (   À8  ]       MPI_GROUP+MPI_CONSTANTS 0   9  H   a   MPI_GROUP%MPI_VAL+MPI_CONSTANTS -   e9  ]       MPI_ERRHANDLER+MPI_CONSTANTS 5   Â9  H   a   MPI_ERRHANDLER%MPI_VAL+MPI_CONSTANTS *   
:  ]       MPI_MESSAGE+MPI_CONSTANTS 2   g:  H   a   MPI_MESSAGE%MPI_VAL+MPI_CONSTANTS '   ¯:  ]       MPI_FILE+MPI_CONSTANTS /   ;  H   a   MPI_FILE%MPI_VAL+MPI_CONSTANTS &   T;  ]       MPI_WIN+MPI_CONSTANTS .   ±;  H   a   MPI_WIN%MPI_VAL+MPI_CONSTANTS '   ù;  ]       MPI_COMM+MPI_CONSTANTS /   V<  H   a   MPI_COMM%MPI_VAL+MPI_CONSTANTS *   <  ]       MPI_REQUEST+MPI_CONSTANTS 2   û<  H   a   MPI_REQUEST%MPI_VAL+MPI_CONSTANTS '   C=  ]       MPI_INFO+MPI_CONSTANTS /    =  H   a   MPI_INFO%MPI_VAL+MPI_CONSTANTS 