#ifndef __STAMPER__H__
#define __STAMPER__H__

#include <vector>

#include <vector>
#include <iostream>
#include <fstream>

#include <cmath>
#include <ASTex/easy_io.h>
#include <Eigen/Core>

#include <ASTex/Stamping/stamp.h>
#include <ASTex/Stamping/sampler.h>

//TODO: consider removing the const of std::vector<const StampBase<I> *>
//if you think it should be authorized to modify the stamps given by the caller.

namespace ASTex
{

namespace Stamping
{


template<typename I>
/**
 * @brief The StamperBase class is an interface for every Stampers.
 * The main idea is to combine the work of a Sampler, one or several Stamps
 * and a Stamp-chosing function to generate texture by stamping the stamps and
 * mixing them in different ways.
 */
class StamperBase
{
public:

	/**
	 * @brief StamperBase constructor of StamperBase.
	 * @param sampler the sampler generating an array of 2D vertices between 0 and 1.
	 * @param stamp a first stamp to add, doesn't add anything if 0.
	 */
	StamperBase(SamplerBase *sampler=0, const StampBase<I> *stamp=0);

	/**
	 * @brief ~StamperBase destructor of StamperBase.
	 * Explicitely desallocates the index generator of the class.
	 */
	virtual ~StamperBase()                                      {delete(m_stampIndexGenerator);}

	//nested types

	typedef typename StampBase<I>::PixelType PixelType;
	/**
	 * @brief The Func_stampIndexGeneratorBase class is a base for an index generator.
	 * An index generator is a functor which choses a stamp index everytime it is called.
	 */
	class Func_stampIndexGeneratorBase
	{
	public:
		Func_stampIndexGeneratorBase() {}
		virtual ~Func_stampIndexGeneratorBase() {}
		virtual size_t operator() (void)=0;
	};

	/**
	 * @brief The Func_stampIndexGeneratorZero class is an override for Func_stampIndexGeneratorBase.
	 * It always returns 0 when called.
	 */
	class Func_stampIndexGeneratorZero : public Func_stampIndexGeneratorBase
	{
	public:
		Func_stampIndexGeneratorZero() : Func_stampIndexGeneratorBase() {}
		size_t operator() (void) {return 0;}
	};

	//main function

	/**
	 * @brief generate a virtual pure method. It is used to generate an output image of size (imageWidth, imageHeight).
	 * @param imageWidth the width of the output image.
	 * @param imageHeight the height of the output image.
	 * @return an Image of the template type I.
	 */
	virtual I generate(int imageWidth, int imageHeight) const = 0;

	//get

	/**
	 * @return the sampler class is used
	 * to generate an array of coordinates during the convolution process.
	 */
	SamplerBase* sampler() const  {return m_sampler;}

	/**
	 * @param index the position in the vector of stamps.
	 * @return a reference on a stamp pointer which may be used during the convolution process.
	 */
	const StampBase<I>* stamp(size_t index)                     {return m_stampVec[index];}

	/**
	 * @brief same as stamp(index), const version (result reference cannot be modified).
	 */
	const StampBase<I>* const & stamp(size_t index) const       {return m_stampVec[index];}

	/**
	 * @brief operator [] same as stamp(index).
	 */
	const StampBase<I>*& operator[] (size_t index)              {return m_stampVec[index];}
	/**
	 * @brief operator [] calls stamp(index), const version (result reference cannot be modified).
	 */
	const StampBase<I>* const & operator[] (size_t index) const {return m_stampVec[index];}

	size_t stampVecSize() const                                 {return m_stampVec.size();}

	//set

	/**
	 * @param sampler is the sampler which will be used
	 * to generate an array of coordinates during the convolution process.
	 */
	void setSampler(const SamplerBase* sampler)                 {m_sampler = sampler;}
	/**
	 * @brief setStamp changes the stamp at vector position index to stamp.
	 * @param stamp this stamp may be used during the convolution process.
	 * @param index the position of the vector on which the stamp is modified.
	 * @pre index is in the range of available stamps.
	 */
	void setStamp(const StampBase<I>* stamp, size_t index)      {m_stampVec[index] = stamp;}

	//add/remove

	/**
	 * @brief addStamp adds a stamp to the collection of stamps.
	 * It is pushed at the end of the internal vector.
	 * @param stamp the stamp to be added.
	 */
	void addStamp(const StampBase<I>* stamp)                    {m_stampVec.push_back(stamp);}
	/**
	 * @brief removeStamp removes a stamp at position index.
	 * Warning: complexity is linear.
	 * @param index position of the stamp.
	 */
	void removeStamp(size_t index)                              {m_stampVec.erase(m_stampVec.begin() + index); }

protected:

	/**
	 * @brief assert_null_before_generate_ crashes the program if the sampler or some stamps are NULL.
	 */
	void assert_null_before_generate_() const;

	SamplerBase                         *m_sampler;
	std::vector<const StampBase<I> *>   m_stampVec;
	Func_stampIndexGeneratorBase        *m_stampIndexGenerator;
};

template<typename I>
/**
 * @brief The StamperBombing class is an override of the StamperBase, used to produce bombing synthesis outputs.
 * The mixing method is to override every written pixels.
 */
class StamperBombing : public StamperBase<I>
{
public:

	/**
	 * @brief StamperBombing constructor of StamperBombing.
	 * @param sampler the sampler generating an array of 2D vertices between 0 and 1.
	 * @param stamp a first stamp to add, doesn't add anything if 0.
	 */
	StamperBombing(SamplerBase *sampler=0, const StampBase<I> *stamp=0);

	/**
	 * @brief generate produces an output image which is the result of a bombing process of a single stamp.
	 * @param imageWidth the width of the output
	 * @param imageHeight the height of the output
	 * @return the output image
	 */
	I generate(int imageWidth, int imageHeight) const;
};

template<typename I>
class StamperTexton : public StamperBase<I>
{
public:

	//warning: input must be a shifted and normalized texton (mean 0).
	/**
	 * @brief StamperTexton constructor of StamperTexton.
	 * @param sampler the sampler generating an array of 2D vertices between 0 and 1.
	 * @param stamp a first stamp to add, doesn't add anything if 0.
	 */
	StamperTexton(SamplerBase *sampler=0, const StampBase<I> *stamp=0);

	//warning: output is a shifted and normalized texton as well.
	/**
	 * @brief generate produces an output image which is the result of a spot noise of a single stamp.
	 * @param imageWidth the width of the output
	 * @param imageHeight the height of the output
	 * @return the output image
	 */
	I generate(int imageWidth, int imageHeight) const;

	//get

	/**
	 * @return whether the output image is periodic or not.
	 * The texton/spot noise is easier to produce when it is a periodic process,
	 * but there exists a workaround with "margins" if it shouldn't be periodic.
	 */
	bool periodicity() {return m_periodicity;}
	/**
	* @return whether margins are used or not.
	* Margins are made to extend the spot noise domain,
	* such that a non periodic image would still be able to receive the full energy
	* of the stamp during the spot noise process.
	* It is strongly advised you let this on to get correct results.
	*/
	bool useMargins() {return m_useMargins;}


	//set

	/**
	 * @param periodicity defines whether the output image should be periodic or not.
	 * The texton/spot noise is easier to produce when it is a periodic process,
	 * but there exists a workaround with "margins" if it shouldn't be periodic.
	 */
	void setPeriodicity(bool periodicity) {m_periodicity = periodicity;}
	/**
	 * @param use defines whether margins should be used or not.
	 * Margins are made to extend the spot noise domain,
	 * such that a non periodic image would still be able to receive the full energy
	 * of the stamp during the spot noise process.
	 * It is strongly advised you let this on to get correct results.
	 */
	void setUseMargins(bool use) {m_useMargins = use;}

private:

	bool m_periodicity; //< if the stamping process is periodic or not.

	bool m_useMargins; //< if true, extends the domain to preserve energy.
};

template<typename I>
StamperBase<I>::StamperBase(SamplerBase *sampler, const StampBase<I> *stamp)
	: m_sampler(sampler), m_stampVec(), m_stampIndexGenerator(0)
{
	if(stamp!=0) //allows the user to insert a default stamp
		m_stampVec.push_back(stamp);
}

template<typename I>
void StamperBase<I>::assert_null_before_generate_() const
{
	assert(m_sampler != NULL && "StamperBase::generate(w, h): sampler uninitialized (use setSampler(s) to give one)");
	for(typename std::vector<const StampBase<I> *>::const_iterator it=m_stampVec.begin(); it!=m_stampVec.end(); ++it)
		assert( (*it) != NULL && "StamperBase::generate(w, h): stamp uninitialized found (use addStamp(s) to give one)");
	return;
}

template<typename I>
StamperBombing<I>::StamperBombing(SamplerBase *sampler, const StampBase<I> *stamp) :
	StamperBase<I>(sampler, stamp)
{
	/// v it doesn't have to look this complicated if you define your own functor class for index generation.
	/// note that it will automatically be deleted by the base class destructor.
	this->m_stampIndexGenerator = new typename StamperBase<I>::Func_stampIndexGeneratorZero;
}

template<typename I>
I StamperBombing<I>::generate(int imageWidth, int imageHeight) const
{
	I im_out;
	int i, j, stampWidth, stampHeight, tx, ty, rx, ry;
	std::vector<Eigen::Vector2f> verticesArray; //stores the result of the sampler
	const StampBase<I> *stamp=0;
	this->assert_null_before_generate_(); /// < you can use this (or not) to check for a null sampler or a null stamp
	verticesArray = this->m_sampler->generate();

	im_out.initItk(imageWidth, imageHeight, true);

	/// v example of iteration over the sampler's result
	for(std::vector<Eigen::Vector2f>::const_iterator it=verticesArray.begin(); it!=verticesArray.end(); ++it)
	{
		stamp = this->m_stampVec[(*this->m_stampIndexGenerator)()]; /// < example of how to utilize a generator

		stampWidth = stamp->width();
		stampHeight = stamp->height();

		i = im_out.width() * (*it)[0]; //i & j: single point coordinates in im_out
		j = im_out.height() * (*it)[1];

		rx = i-stampWidth/2; //region origin coordinates (int)
		ry = j-stampHeight/2;
		Region reg = gen_region(rx, ry, stampWidth, stampHeight); /// < example of how to utilize regions for stamping
		im_out.for_region_pixels(reg, [&] (typename StampBase<I>::PixelType& pix, int x, int y)
		{
			if(x>0 && x<im_out.width() && y>0 && y<im_out.height())
			{
				tx=x-rx; //stamp coordinate shifted from the region origin coordinates
				ty=y-ry;
				pix = stamp->pixel(tx, ty); /// < example of mixing function (here, erase previous pixel with new one)
			}
		});
	}

	return im_out;
}

template<typename I>
StamperTexton<I>::StamperTexton(SamplerBase *sampler, const StampBase<I> *stamp) :
	StamperBase<I>(sampler, stamp),
	m_periodicity(false),
	m_useMargins(true)
{
	this->m_stampIndexGenerator = new typename StamperBase<I>::Func_stampIndexGeneratorZero;
}

template<typename I>
I StamperTexton<I>::generate(int imageWidth, int imageHeight) const
{
	I im_out;
	std::vector<Eigen::Vector2f> verticesArray;
	double stampWidth=0.0, stampHeight=0.0;
	double i, j;
	double otx, oty; //texton origin in texture space (top left)
	double tx, ty; //texton coordinates
	double nbHit, nbHitPerPixel; //keeps track of the number of times the texture was hit
	double lambda; //lambda parameter in poisson process
	//assuming given texton is a zero mean texture with normalized variance

	const StampBase<I> *stamp=0;

	this->assert_null_before_generate_();
	verticesArray = this->m_sampler->generate();

	im_out.initItk(imageWidth, imageHeight, true);

	nbHit=0;

	for(std::vector<Eigen::Vector2f>::const_iterator it=verticesArray.begin(); it!=verticesArray.end(); ++it)
	{
		stamp = this->m_stampVec[(*this->m_stampIndexGenerator)()];
		stampWidth = stamp->width();
		stampHeight = stamp->height();

		if(m_periodicity || (!m_periodicity && !m_useMargins))
		{
			i = imageWidth * (*it)[0];
			j = imageHeight * (*it)[1];
		}
		else
		{
			//we need to introduce margins here, to shoot textons outside of the domain
			i = (imageWidth + stampWidth ) * (*it)[0] - stampWidth/2.0;
			j = (imageHeight + stampHeight ) * (*it)[1] - stampHeight/2.0;
		}

		if(m_periodicity)
		{
			otx=i; //texton origin in texture space (top left)
			oty=j;
		}
		else
		{
			otx=i-stampWidth/2.0; //texton origin in texture space (top left)
			oty=j-stampHeight/2.0;
		}

		Region reg = gen_region(std::floor(otx), std::floor(oty), stampWidth+1, stampHeight+1); //note: regions are weak when shooting between pixels

		nbHit += stampWidth*stampHeight; //with periodicity, the entire energy of the texton hits the texture all the time.
		//Without, we pretend it did and normalize including the size of the margins.
		if(m_periodicity)
		{
			im_out.for_region_pixels(reg, [&] (typename I::PixelType& pix, int x, int y) //with periodicity
			{
				(void)pix; /// pix is unused: remove parameter if there can be a (int, int) lambda
				tx=x-otx; //texton coordinate in texton space
				ty=y-oty; //texton coordinate

				im_out.pixelAbsolute((x+im_out.width())%imageWidth, (y+im_out.height())%imageHeight) += stamp->pixel(tx, ty);
			});
		}
		else
		{

			//the region we stamp : it is one pixel longer (per dim.) when the texton can stamped between pixels <=> when bilinear interpolation is activated
			//in othger terms, the region is a 2D bounding box for the texton
			im_out.for_region_pixels(reg, [&] (typename I::PixelType& pix, int x, int y) //without periodicity
			{
				(void)pix; /// pix is unused: remove parameter if there can be a (int, int) lambda
				if(x>=0 && y>=0 && x<imageWidth && y<imageHeight)
				{ //here we compute pixel values, it's additive given there is 0 outside of the tx range.
					tx=x-otx; //texton coordinate in texton space
					ty=y-oty; //texton coordinate

					im_out.pixelAbsolute(x, y) += stamp->pixel(tx, ty);

				}
			});
		}
	}
	//warning: assumes we use only one stamp / stamps with same width and height
	nbHitPerPixel = m_periodicity || !m_useMargins ? nbHit/(imageWidth*imageHeight) : nbHit/((imageWidth+stampWidth)*(imageHeight+stampHeight));

	lambda = double(nbHitPerPixel)/(stampWidth * stampHeight);

	std::cout << "Texton noise: number of impacts per pixel: " << std::to_string(nbHitPerPixel) << std::endl;

	im_out.for_all_pixels([&] (typename I::PixelType &pix) {
			pix = pix * (1.0/sqrt(lambda)) /*+ mean*/;
	});

	return im_out;
}

} //namespace Stamping

} //namespace ASTex


#endif //__STAMPER__H__
