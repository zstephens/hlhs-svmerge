
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

def subtractive_filter(myList, items_to_remove):
	outList = []
	j0 = 0
	for i in xrange(len(myList)):
		toRemove = False
		while items_to_remove[j0][0][1] <= myList[i][0][0]:
			j0 += 1
		for j in xrange(j0,len(items_to_remove)):
			if could_be_same_sv(myList[i][1],items_to_remove[j][1]):
				toRemove = True
				break
			if items_to_remove[j][0][0] > myList[i][0][1]:
				break
		if toRemove == False:
			outList.append(myList[i])
	return outList


def intersection_filter(myList1, myList2):
	outList = []
	j0 = 0
	for i in xrange(len(myList)):
		isFound = False
		while myList2[j0][0][1] <= myList[i][0][0]:
			j0 += 1
		for j in xrange(j0,len(myList2)):
			if could_be_same_sv(myList[i][1],myList2[j][1]):
				isFound = True
				break
			if myList2[j][0][0] > myList[i][0][1]:
				break
		if isFound == True:
			outList.append(myList[i])
	return outList


