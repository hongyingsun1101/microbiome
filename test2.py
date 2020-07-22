def plot_roc_curve_microbiome_test(pplacer_ref_list, pplacer_stats_list, community_list, cutoff_list, test_data_list, scoreOption=True, testOption=False):
    for refIndex in range(len(pplacer_ref_list)):
        pplacer_ref = pplacer_ref_list[refIndex]
        for statsIndex in range(len(pplacer_stats_list)):
            pplacer_stats = pplacer_stats_list[statsIndex]
            for communityIndex in range(len(community_list)):
                community = community_list[communityIndex]
                for i in range(len(cutoff_list)):
                    cutoff=cutoff_list[i]
                    if(is_float(cutoff)):
                        cutoff_binary=float(cutoff)
                    else:
                        if(scoreOption):
                            cutoff_binary=float(df[pplacer_ref+community].describe().loc[[cutoff]])

                        else:
                            cutoff_binary = float(df[pplacer_ref+pplacer_stats].describe().loc[[cutoff]])
                    # no test situation, which is the default option
                    if (not testOption):

                        if(scoreOption):
                            mask = df[pplacer_ref+community] <=  cutoff_binary
                            df.loc[mask, pplacer_ref+community+'_binary'] = 1
                            mask = df[pplacer_ref+community] >cutoff_binary
                            df.loc[mask, pplacer_ref+community+'_binary'] = 0
                            df_binary = df[[pplacer_ref+pplacer_stats, pplacer_ref+community+'_binary']].dropna()
                            data_stats = df_binary[pplacer_ref+pplacer_stats].to_numpy().reshape(-1,1)
                            binary_label =  df_binary[pplacer_ref+community+'_binary'].to_numpy()
                            print(' The score cutoff '+ cutoff +' for Reference ' + pplacer_ref +' community ' + community   + ' with pplacer_stats '+ pplacer_stats[1:] + ': %.2f' % cutoff_binary )
                            # plot_roc(data_stats,binary_label)
                            plot_roc_microbiome(data_stats,binary_label,x=None,y=None,data_test=False)
                        else:
                            mask = df[pplacer_ref+pplacer_stats] <=  cutoff_binary
                            df.loc[mask, pplacer_ref+pplacer_stats+'_binary'] = 1
                            mask = df[pplacer_ref+pplacer_stats] >cutoff_binary
                            df.loc[mask, pplacer_ref+pplacer_stats+'_binary'] = 0
                            df_binary = df[[pplacer_ref+community, pplacer_ref+pplacer_stats+'_binary']].dropna()
                            data_stats = df_binary[pplacer_ref+community].to_numpy().reshape(-1,1)
                            binary_label =  df_binary[pplacer_ref+pplacer_stats+'_binary'].to_numpy()
                            print(' The pplacer_stats_cutoff '+ cutoff +' for Reference ' + pplacer_ref +' community ' + community + ' pplacer_stats '  + pplacer_stats[1:]  + ': %.2f' % cutoff_binary )
                            # plot_roc(data_stats,binary_label)
                            plot_roc_microbiome(data_stats,binary_label,x=None,y=None,data_test=False)
                    
                    # the no test situation
                    else:
                        for j in range(len(test_data_list)):
                            test=test_data_list[j]
                            if(scoreOption):
                                mask = df[pplacer_ref+community] <=  cutoff_binary
                                df.loc[mask, pplacer_ref+community+'_binary'] = 1
                                mask = df[pplacer_ref+community] >cutoff_binary
                                df.loc[mask, pplacer_ref+community+'_binary'] = 0
                                df_binary = df[[pplacer_ref+pplacer_stats, pplacer_ref+community+'_binary']].dropna()
                                data_stats = df_binary[pplacer_ref+pplacer_stats].to_numpy().reshape(-1,1)
                                binary_label =  df_binary[pplacer_ref+community+'_binary'].to_numpy()
                                
                                mask_test = df[test+community] <=  cutoff_binary
                                df.loc[mask_test, test+community+'_binary'] = 1
                                mask_test = df[test+community] >cutoff_binary
                                df.loc[mask_test, test+community+'_binary'] = 0
                                df_binary = df[[test+pplacer_stats, test+community+'_binary']].dropna()
                                x = df_binary[test+pplacer_stats].to_numpy().reshape(-1,1)
                                y =  df_binary[test+community+'_binary'].to_numpy()

                                print(' The score cutoff '+ cutoff +' for Reference ' + pplacer_ref +' community ' + community   + ' with pplacer_stats '+ pplacer_stats[1:] + 'compared with test ' + test +': %.2f' % cutoff_binary )
                                plot_roc_microbiome(data_stats,binary_label,x,y,data_test=True)
                            else:

                                mask = df[pplacer_ref+pplacer_stats] <=  cutoff_binary
                                df.loc[mask, pplacer_ref+pplacer_stats+'_binary'] = 1
                                mask = df[pplacer_ref+pplacer_stats] >cutoff_binary
                                df.loc[mask, pplacer_ref+pplacer_stats+'_binary'] = 0
                                df_binary = df[[pplacer_ref+pplacer_stats, pplacer_ref+pplacer_stats+'_binary']].dropna()
                                data_stats = df_binary[pplacer_ref+community].to_numpy().reshape(-1,1)
                                binary_label =  df_binary[pplacer_ref+pplacer_stats+'_binary'].to_numpy()
                                
                                mask_test = df[test+pplacer_stats] <=  cutoff_binary
                                df.loc[mask_test, test+pplacer_stats+'_binary'] = 1
                                mask_test = df[test+pplacer_stats] >cutoff_binary
                                df.loc[mask_test, test+pplacer_stats+'_binary'] = 0
                                df_binary = df[[test+pplacer_stats, test+pplacer_stats+'_binary']].dropna()
                                x = df_binary[test+community].to_numpy().reshape(-1,1)
                                y =  df_binary[test+pplacer_stats+'_binary'].to_numpy()
                                binary_label =  df_binary[pplacer_ref+pplacer_stats+'_binary'].to_numpy()
                                print(' The pplacer_stats_cutoff '+ cutoff +' for Reference ' + pplacer_ref +' community ' + community + ' pplacer_stats '  + pplacer_stats[1:] + 'compared with test ' + test  + ': %.2f' % cutoff_binary )
                                plot_roc_microbiome(data_stats,binary_label,x,y,data_test=True)