def plot_roc_curve_microbiome_2(pplacer_ref_list, pplacer_stats_list, community_list, cutoff_list, scoreOption=True):
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
                    if(scoreOption):
                        mask = df[pplacer_ref+community] <=  cutoff_binary
                        df.loc[mask, pplacer_ref+community+'_binary'] = 1
                        mask = df[pplacer_ref+community] >cutoff_binary
                        df.loc[mask, pplacer_ref+community+'_binary'] = 0
                        df_binary = df[[pplacer_ref+pplacer_stats, pplacer_ref+community+'_binary']].dropna()
                        data_stats = df_binary[pplacer_ref+pplacer_stats].to_numpy().reshape(-1,1)
                        binary_label =  df_binary[pplacer_ref+community+'_binary'].to_numpy()
                        print(' The score cutoff '+ cutoff +' for Reference ' + pplacer_ref +' community ' + community   + ' with pplacer_stats '+ pplacer_stats[1:] + ': %.2f' % cutoff_binary )
                        plot_roc(data_stats,binary_label)
                    else:
                        mask = df[pplacer_ref+pplacer_stats] <=  cutoff_binary
                        df.loc[mask, pplacer_ref+pplacer_stats+'_binary'] = 1
                        mask = df[pplacer_ref+pplacer_stats] >cutoff_binary
                        df.loc[mask, pplacer_ref+pplacer_stats+'_binary'] = 0
                        df_binary = df[[pplacer_ref+community, pplacer_ref+pplacer_stats+'_binary']].dropna()
                        data_stats = df_binary[pplacer_ref+community].to_numpy().reshape(-1,1)
                        binary_label =  df_binary[pplacer_ref+pplacer_stats+'_binary'].to_numpy()
                        print(' The pplacer_stats_cutoff '+ cutoff +' for Reference ' + pplacer_ref +' community ' + community + ' pplacer_stats '  + pplacer_stats[1:]  + ': %.2f' % cutoff_binary )
                        plot_roc(data_stats,binary_label)
                


def plot_roc_curve(fpr, tpr):
    plt.plot(fpr, tpr, color='orange', label='ROC')
    plt.plot([0, 1], [0, 1], color='darkblue', linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.legend()
    plt.show()


def plot_roc(data_X, class_label):
    trainX, testX, trainy, testy = train_test_split(data_X, class_label, test_size=0.3, random_state=1)
    model = RandomForestClassifier()
    model.fit(trainX, trainy)
    probs = model.predict_proba(testX)
    probs = probs[:, 1]
    auc = roc_auc_score(testy, probs)
    fpr, tpr, thresholds = roc_curve(testy, probs)
    print('AUC: %.2f' % auc)
#     print( thresholds)
#     print('Model: ')
#     print(model)
    plot_roc_curve(fpr, tpr)



