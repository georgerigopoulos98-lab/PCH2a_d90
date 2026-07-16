# Voxhunt Allen refferences:
# developing mouse E18.5, P4 and P14
# developing human 21 pcw cerebrum and brainstem

#Before exporting, verify structure:
print(adata)
print(adata.obs_names[:5]) #They should look like: Index(['AAACCTGAG...', 'AAACCTGAG...', ...]) and not: Index([0,1,2,3,...]).
print(adata.var_names[:5])
