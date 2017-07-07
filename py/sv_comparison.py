
# how far away can SV call positions be such that we'll still consider them the same event?
ENDPOINT_BUFFER = 150

def could_be_same_sv(sv_cluster1, sv_cluster2):
	# is same event type
	if sv_cluster1[0][0] != sv_cluster2[0][0]:
		return False

	# on same chr
	if sv_cluster1[0][1] != sv_cluster2[0][1]:
		return False

	# some start/end positions within each cluster are close enough together
	valid_start = []
	for i in xrange(len(sv_cluster1)):
		for j in xrange(len(sv_cluster2)):
			if abs(sv_cluster1[i][2] - sv_cluster2[j][2]) <= ENDPOINT_BUFFER:
				if abs(sv_cluster1[i][3] - sv_cluster2[j][3]) <= ENDPOINT_BUFFER:
					valid_start.append((i,j))
	if len(valid_start) == 0:
		return False

	# congrats, you passed!
	return True

#
#	ASSUMES INPUT LISTS ARE SORTED
#
def sv_filter(list1,list2,op='subtract'):
	if op not in ['subtract','intersect']:
		print '\nError: filter op must be "subtract" or "intersect"\n'
		exit(1)
	outList = []
	chrMap  = {}
	for i in xrange(len(list1)):
		mc1 = list1[i][0][0]
		if mc1 not in chrMap:
			chrMap[mc1] = [i,None]
	for i in xrange(len(list2)):
		mc2 = list2[i][0][0]
		if mc2 in chrMap and chrMap[mc2][1] == None:
			chrMap[mc2][1] = i

	for i in xrange(len(list1)):
		j0    = chrMap[list1[i][0][0]][1]
		in_l2 = False
		if j0 != None:
			for j in xrange(j0,len(list2)):
				if list2[j][0][0] != list1[i][0][0] or list2[j][0][1] > list1[i][0][2]:
					break
				if could_be_same_sv(list1[i][1],list2[j][1]):
					in_l2 = True
					break
		if op == 'intersect' and in_l2:
			outList.append(list1[i])
		elif op == 'subtract' and not(in_l2):
			outList.append(list1[i])

	return outList


