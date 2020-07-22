fig, ax1=plt.subplots()
ax1.hist([df.A_edpl, df.B_edpl,df.C_edpl, df.D_edpl,df.E_edpl], color=['b','g','r','c','m'],label=['A','B','C','D','E'])
ax1.set_ylim(0,4000)
ax1.set_ylabel("Count")
plt.title("edpl")
plt.legend(loc='upper right')
plt.tight_layout()
# plt.show()
plt.savefig('edpl.png',dpi=1000)



ax1.set_xlim(0,max(df.A_prichness.max(), df.B_prichness.max(),df.C_prichness.max(), df.D_prichness.max(), df.E_prichness.max()))


stats.percentileofscore(df.A_prichness,5),stats.percentileofscore(df.B_prichness,5),stats.percentileofscore(df.C_prichness,5),stats.percentileofscore(df.D_prichness,5),stats.percentileofscore(df.E_prichness,5)