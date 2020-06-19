#ifndef CHROMAP_CACHE_H_
#define CHROMAP_CACHE_H_

#include "index.h"


namespace chromap {

struct _mm_cache_entry
{
	std::vector<uint64_t> minimizers ;
	std::vector<int> offsets ; // the distance to the next minimizer
	std::vector<struct _candidate> positive_candidates ;
	std::vector<struct _candidate> negative_candidates ;

	int weight ;
} ;

class mm_cache
{
private:
	int cache_size ;
	struct _mm_cache_entry *cache ;
	
	// 0: not match. -1: opposite order. 1: same order
	int IsMinimizersMatchCache(const std::vector<std::pair<uint64_t, uint64_t> > &minimizers, const struct _mm_cache_entry &cache)
	{
		if (cache.minimizers.size() != minimizers.size())
			return 0 ;
		int size = minimizers.size() ;
		int i, j ;
		int direction = 0 ;
		for (i = 0 ; i < size ; ++i)
		{
			if (cache.minimizers[i] != minimizers[i].first)
				break ;
		}
		if (i >= size)
		{
			for (i = 0 ; i < size - 1 ; ++i)	
			{
				if (cache.offsets[i] != ((int)(minimizers[i + 1].second)>>1) - ((int)(minimizers[i].second)>>1))
					break ;
			}
			if (i >= size - 1)
				direction = 1 ;
		}
		
		if (direction == 1)
			return 1 ;
		
		for (i = 0, j = size - 1; i < size; ++i, --j)
		{
			if (cache.minimizers[i] != minimizers[j].first)
				break ;
		}
		if (i >= size)
		{
			for (i = 0, j = size - 1; i < size - 1; ++i, --j)
			{
				if (cache.offsets[i] != ((int)(minimizers[j].second)>>1) - ((int)(minimizers[j - 1].second)>>1))
					break ;
			}

			if (i >= size - 1)
			{
				direction = -1 ;
			}
		}

		return direction ;
	}
public:
	mm_cache(int size) 
	{
		cache = new struct _mm_cache_entry[size] ;
		cache_size = size ;
		memset(cache, 0, sizeof(cache[0]) * size) ;
	}
	~mm_cache()
	{
		delete[] cache ;
	}
	
	// Return the hash entry index. -1 if failed.
	int Query(const std::vector<std::pair<uint64_t, uint64_t> > &minimizers, 
			std::vector<struct _candidate> &pos_candidates, std::vector<struct _candidate> &neg_candidates, 
			uint32_t read_len)
	{
		int i ;
		int msize = minimizers.size() ;
		uint64_t h = 0 ;
		for (i = 0 ; i < msize; ++i)
			h ^= (minimizers[i].first >> 8) ;	
		int hidx = h % cache_size ;
		int direction = IsMinimizersMatchCache(minimizers, cache[hidx]) ;
		if (direction == 1)
		{
			pos_candidates = cache[hidx].positive_candidates ;
			neg_candidates = cache[hidx].negative_candidates ;
			return hidx ;
		}
		else if (direction == -1)
		{
			int size = cache[hidx].negative_candidates.size() ;
			pos_candidates = cache[hidx].negative_candidates ;
			for (i = 0 ; i < size ; ++i)
				pos_candidates[i].refPos = cache[hidx].negative_candidates[i].refPos - read_len + 1 ;

			size = cache[hidx].positive_candidates.size() ;
			neg_candidates = cache[hidx].positive_candidates ;
			for (i = 0 ; i < size ; ++i)
				neg_candidates[i].refPos = cache[hidx].positive_candidates[i].refPos + read_len - 1 ;
			return hidx ;
		}
		else
		{
			return -1 ;
		}
	}

	void Update(const std::vector<std::pair<uint64_t, uint64_t> > &minimizers, 
			std::vector<struct _candidate> &pos_candidates, std::vector<struct _candidate> &neg_candidates )
	{
		int i ;
		int msize = minimizers.size() ;

		uint64_t h = 0 ;
		for (i = 0 ; i < msize; ++i)
			h ^= (minimizers[i].first >> 8) ;	
		int hidx = h % cache_size ;

		int direction = IsMinimizersMatchCache(minimizers, cache[hidx]) ;
		if (direction != 0)
			++cache[hidx].weight ;
		else
			--cache[hidx].weight ;
		// Renew the cache
		if (cache[hidx].weight <= 0 )
		{
			cache[hidx].weight = 1 ;
			int size = minimizers.size() ;
			cache[hidx].minimizers.resize(size) ;
			cache[hidx].offsets.resize(size - 1) ;
			for (i = 0 ; i < size ; ++i)
				cache[hidx].minimizers[i] = minimizers[i].first ;
			for (i = 0 ; i < size - 1; ++i)
			{
				cache[hidx].offsets[i] = ((int)(minimizers[i + 1].second)>>1) - ((int)(minimizers[i].second)>>1) ;
			}
			cache[hidx].positive_candidates = pos_candidates ;
			cache[hidx].negative_candidates = neg_candidates ;
		}
	}
	
		
} ;

} 

#endif
