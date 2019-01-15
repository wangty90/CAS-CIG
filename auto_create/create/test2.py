#!/usr/bin/python
mrg_init_ccc=[['atm','a','a','a',['l','o','i'],'xao_ax'],\
              ['wrf','w','m','m',['c'],2],\
              ['cam','c','c','c',['m']],\
              ['geam','ge','ge','chem',['ca'],2],\
              ['cam','ca','ca','ca',['chem']],\
              ['ice','i','i','i',['a','o']],\
              ['ocn','o','o','o',['a','i','r']],\
              ['lnd','l','l','l',['a']],\
              ['glc','g','g','g',['s']],\
              ['sno','s','s','s',['g']]\
]
clist = mrg_init_ccc
for c in range(len(clist)):
    print clist[c][4],len(clist[c])
    for cn in range(len(clist[c][4])):
	print clist[c][4][cn]

