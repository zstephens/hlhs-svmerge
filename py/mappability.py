import bisect

# mappability track
class MappabilityTrack:
	def __init__(self, bedfile):
		#print 'reading input mappability track...'
		f = open(bedfile,'r')
		self.all_tracks = {}
		for line in f:
			splt = line.strip().split('\t')
			[myChr,myPos,myEnd] = [splt[0],int(splt[1]),int(splt[2])]
			if myChr not in self.all_tracks:
				self.all_tracks[myChr] = [-1]
			self.all_tracks[myChr].append(myPos)
			self.all_tracks[myChr].append(myEnd)
		f.close()

	# return True if index is in bedfile region, Flase if outside
	def query(self, myChr, coord, endPointInclusive=True):
		if myChr in self.all_tracks:
			bcoord = bisect.bisect(self.all_tracks[myChr],coord)
			if bcoord%2 == 0:
				return True
			if endPointInclusive and coord == self.all_tracks[myChr][bcoord-1]:
				return True
		return False

	# returns number of positions of a range that lie in bedfile regions
	def query_range(self, myChr, coord1, coord2, query_endPointInclusive=False):
		count = 0
		if myChr in self.all_tracks:
			for i in xrange(coord1,coord2+1*query_endPointInclusive):
				if self.query(myChr,i):
					count += 1
		return count

	# more efficient version of the above...
	def query_range_faster(self, myChr, coord1, coord2, query_endPointInclusive=False):
		count = 0
		if myChr in self.all_tracks:
			lb = bisect.bisect(self.all_tracks[myChr],coord1)
			ub = bisect.bisect(self.all_tracks[myChr],coord2+1*query_endPointInclusive)
			
			# entirely in a single region
			if lb == ub and lb%2 == 0:
				return coord2 - coord1 + 1*query_endPointInclusive
			if lb == ub and lb%2 == 1:
				return 0

			#print [coord1,coord2],(lb,ub)
			size_list   = [self.all_tracks[myChr][lb]-coord1]
			in_out_list = [(lb%2 == 0)]
			for i in xrange(lb+1,ub):
				size_list.append(self.all_tracks[myChr][i]-self.all_tracks[myChr][i-1])
				in_out_list.append((i%2 == 0))
			size_list.append(coord2-self.all_tracks[myChr][ub-1])
			in_out_list.append((ub%2 == 0))

			for i in xrange(len(size_list)):
				if in_out_list[i]:
					count += size_list[i] + 1*query_endPointInclusive

		return count