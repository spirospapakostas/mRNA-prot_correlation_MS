## 040518 - this scipt is for replicating results in Papakostas et al 2018 related to estimation of QST and thermal variation explained in protein and mRNA data, corresponds to Figure3, Figure S4, Figure S5.

    # PREPERAING DATA
    
        #loading all gene expresison data. these values are the normalised values as descibed in the MS
        #note that the number of gene in this list that corresponses protein IDs is slightly higher than 654, due to the fact that some protein IDs are represneted by multiple mRNA IDs (se Papakostas et al 2018 for details). 
        #note that genes with a protein data expression is indicated in the "X654genes" column. The first part of the ID (before the first "|" indicator provides the ID of the protein that matching IDs betwene proteins and genes)
         genes = read.table ( "file:///C:/Users/Localadmin_aykanat/Desktop/spiros160418/16622genes_w654.txt" , h=T, sep = "\t" , stringsAsFactors = F )
        #loading protein data (note that there are two ID columns. alt.ID is NA id there is only one corresponding mRNA to the protein seqeunce)
         prots = read.table ("C:/Users/Localadmin_aykanat/Desktop/spiros160418/654genes_proteinID.txt" , h=T, sep = "\t" , stringsAsFactors = F )
    
            # removing first two columns and trasposing the data farme
              genes2 = t(genes[,-c(1:3)])
              prots2 = t(prots[,-c(1:3)])
            
            # load files that give each individual thermal, population origin, and the thermal treatment
              load(file = "KEYt")
              load(file = "KEYp")
    
    # Calculating coefficients to estimate QST and % thermal variance explained, sepertately at 10C and 6C thernal treatement
       
        library(lme4)
          ## QST
            
            ## mRNA
             ALL_T_QSTs_10 = sapply ( 1 : ncol(genes2) , function(i) { print(i)
                    RES = genes2[ KEYt$TEMP == 10,i]
                    POP = KEYt$POPS[KEYt$TEMP == 10]
                    A = lmer( RES ~ 1 + (1|POP) )
                    A
                    among = VarCorr(A)$POP[1]
                    within =  attr(VarCorr(A),"sc")^2
                    c(within,among ,among / (among + (2 * within) ))
                    })
             ALL_T_QSTs_6 = sapply ( 1 : ncol(genes2) , function(i) { print(i)
                    RES = genes2[ KEYt$TEMP == 6,i]
                    POP = KEYt$POPS[KEYt$TEMP == 6]
                    A = lmer( RES ~ 1 + (1|POP) )
                    A
                    among = VarCorr(A)$POP[1]
                    within =  attr(VarCorr(A),"sc")^2
                    c(within,among ,among / (among + (2 * within) ))
                    })
               all_Tqst_6 = ALL_T_QSTs_6 [3,] 
               all_Tqst_10 = ALL_T_QSTs_10 [3,] 
            
            ## protein
             P_QSTs_6  = sapply ( 1 : ncol(prots2) , function(i) { print(i)
                    RES = prots2[ KEYp$TEMP == 6,i]
                    POP = KEYp$POPS[KEYp$TEMP == 6]
                    A = lmer( RES ~ 1 + (1|POP) )
                    among = VarCorr(A)$POP[1]
                    within =  attr(VarCorr(A),"sc")^2
                    c(within,among ,among / (among + (2 * within) ))
                    })
             P_QSTs_10 = sapply ( 1 : ncol(prots2) , function(i) { print(i)
                    RES = prots2[ KEYp$TEMP == 10,i]
                    POP = KEYp$POPS[KEYp$TEMP == 10]
                    A = lmer( RES ~ 1 + (1|POP) )
                    among = VarCorr(A)$POP[1]
                    within =  attr(VarCorr(A),"sc")^2
                    c(within,among ,among / (among + (2 * within) ))
                    })
              Pqst_6 = P_QSTs_6 [3,] 
              Pqst_10 = P_QSTs_10 [3,] 
             
        
        
          ## % thermal variance explained
          ## note that thermal origin is confounded completely within populations and therfore correlation of model residuals with and without thermal effect is > 0.999, and all variance accountred for thermal effects can be attributed to population effect.
              
              #mRNA
                 ALLgenes_6T = t(sapply ( 1 : ncol(genes2) , function(i) { print(i)
                      
                      RES = genes2[ KEYt$TEMP == 6,i]
                      POP = KEYt$POPS[KEYt$TEMP == 6]
                      TERMAL = gsub("KV","C",gsub("VA","C",gsub("HA","W",gsub("OT","W",POP))))
            
                      A = lmer( RES ~ 1 + (1|POP)  , REML=T)
                      B = lmer( RES ~ TERMAL + (1|POP)  , REML=T)
            
                      within =  attr(VarCorr(A),"sc")^2
                      withinB =  attr(VarCorr(B),"sc")^2
                      among = VarCorr(A)$POP[1]
                      amongB = VarCorr(B)$POP[1]
                      varTERMAL = among - amongB
                      c(varTERMAL,amongB,withinB,within)
                 }))
                 ALLgenes_10T = t(sapply ( 1 : ncol(genes2) , function(i) { print(i)
                      RES = genes2[ KEYt$TEMP == 10,i]
                      POP = KEYt$POPS[KEYt$TEMP == 10]
                      TERMAL = gsub("KV","C",gsub("VA","C",gsub("HA","W",gsub("OT","W",POP))))
                      
                      A = lmer( RES ~ 1 + (1|POP) , REML=T)
                      B = lmer( RES ~ TERMAL + (1|POP)  , REML=T)
            
                      within =  attr(VarCorr(A),"sc")^2
                      withinB =  attr(VarCorr(B),"sc")^2
                      among = VarCorr(A)$POP[1]
                      amongB = VarCorr(B)$POP[1]
                      varTERMAL = among - amongB
                      c(varTERMAL,amongB,withinB,within)
                 }))
            
                 TOTVAR_10 = apply( ALLgenes_10T[,1:3] ,1 , sum)
                 TOTVAR_6 = apply( ALLgenes_6T[,1:3] ,1 , sum)
                  
                 var_termal_originT6_all  =  ( ALLgenes_6T[,1] / TOTVAR_6 )
                 var_termal_originT10_all =  ( ALLgenes_10T[,1] / TOTVAR_10 )
              
                 ## negative variances equate to zero.
                 var_termal_originT6_all [ var_termal_originT6_all < 0 ] = 0
                 var_termal_originT10_all [ var_termal_originT10_all < 0 ] = 0
            
              #protein
            
                 ALLgenes_10P = t(sapply ( 1 : ncol(prots2) , function(i) { print(i)
                      RES = prots2[ KEYp$TEMP == 10,i]
                      POP = KEYp$POPS[KEYp$TEMP == 10]
                      TERMAL = gsub("KV","C",gsub("VA","C",gsub("HA","W",gsub("OT","W",POP))))
                      
                      A = lmer( RES ~ 1 + (1|POP) , REML=T)
                      B = lmer( RES ~ TERMAL + (1|POP)  , REML=T)
                      C = lmer( RES ~ 1 + (1|TERMAL)  , REML=T)
            
                      within =  attr(VarCorr(A),"sc")^2
                      withinB =  attr(VarCorr(B),"sc")^2
                      among = VarCorr(A)$POP[1]
                      amongB = VarCorr(B)$POP[1]
                      varTERMAL = among - amongB
                      c(varTERMAL,amongB,withinB,within)
                 }))
                 ALLgenes_6P = t(sapply ( 1 : ncol(prots2) , function(i) { print(i)
                      
                      RES = prots2[ KEYp$TEMP == 6,i]
                      POP = KEYp$POPS[KEYp$TEMP == 6]
                      TERMAL = gsub("KV","C",gsub("VA","C",gsub("HA","W",gsub("OT","W",POP))))
            
                      A = lmer( RES ~ 1 + (1|POP)  , REML=T)
                      B = lmer( RES ~ TERMAL + (1|POP)  , REML=T)
                      C = lmer( RES ~ 1 + (1|TERMAL)  , REML=T)
            
                      within =  attr(VarCorr(A),"sc")^2
                      withinB =  attr(VarCorr(B),"sc")^2
                      among = VarCorr(A)$POP[1]
                      amongB = VarCorr(B)$POP[1]
                      varTERMAL = among - amongB
                      c(varTERMAL,amongB,withinB,within)
                 }))
                 
                 TOTVAR_10P = apply( ALLgenes_10P[,1:3] ,1 , sum)
                 TOTVAR_6P = apply( ALLgenes_6P[,1:3] ,1 , sum)
                 
                 var_termal_originP6_all  =  ( ALLgenes_6P[,1] / TOTVAR_6P )
                 var_termal_originP10_all =  ( ALLgenes_10P[,1] / TOTVAR_10P )
                 
                 ## variance less than zero is equelies to zero.
                 var_termal_originP6_all [ var_termal_originP6_all < 0 ] = 0
                 var_termal_originP10_all [ var_termal_originP10_all < 0 ] = 0    
## ENDS         
         
         
