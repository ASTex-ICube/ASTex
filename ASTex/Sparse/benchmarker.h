#ifndef __BENCHMARKER__H__
#define __BENCHMARKER__H__

#include "dictionaryProcessor.h"
#include "utils.h"

/////////////////////

class SparseBenchmarker : public ASTex::Benchmarker
{
public:
	SparseBenchmarker();

	using ImageType = ASTex::ImageRGBd;

	class CompareRGBSum
	{
	public:
		CompareRGBSum() {}
		bool operator()(const ASTex::ImageRGBd::PixelType &pix, const ASTex::ImageRGBd::PixelType &other)
		{
			double normpix = pix[0] + pix[1] + pix[2];
			double normother = other[0] + other[1] + other[2];
			return normpix < normother;
		}
	};

	//step 1: requirements//

	void setInput(const ImageType &input) {m_input=input;}
	void setInputName(const std::string &name) {m_inputName = name;}
	void setLearningDirectory(const std::string &learningDirectory);
	DictionaryProcessor<ImageType> &dictionaryProcessor() {return m_dictionaryProcessor;}

	//step 2: options//

	void setNbResults(unsigned results) {m_numberResults = results;}
	void setNumberIterationsLearning(unsigned nbIterations) {m_numberIterationsLearning = nbIterations;}
	void setNumberIterationsSynthesising(unsigned nbIterations) {m_numberIterationsSynthesising = nbIterations;}
	void setLoadDictionaryMode(bool b) {m_loadDictionaryMode = b;}

	//step 3

	void generate();

private:

	bool			m_loadDictionaryMode;

	ImageType		m_input;

	unsigned m_numberIterationsLearning;
	unsigned m_numberIterationsSynthesising;
	unsigned m_numberResults;

	//result directories

	std::string m_learningDirectory;
	std::string m_inputName;

	DictionaryProcessor<ImageType> m_dictionaryProcessor;
};

SparseBenchmarker::SparseBenchmarker() :
	Benchmarker(),
	m_loadDictionaryMode(false),
	m_input(),
	m_numberIterationsLearning(8),
	m_numberIterationsSynthesising(8),
	m_numberResults(1),
	m_learningDirectory(),
	m_inputName(),
	m_dictionaryProcessor()
{}

void SparseBenchmarker::setLearningDirectory(const std::string &learningDirectory)
{
	m_learningDirectory = learningDirectory;
	ASTex::create_directory(m_root + learningDirectory);
}

void SparseBenchmarker::generate()
{
	std::string ioPath = m_root + m_learningDirectory + "/" + m_inputName
			+ "_learn" + std::to_string(m_numberIterationsLearning)
			+ "_" + std::to_string(m_dictionaryProcessor.nbAtoms()) + "atoms"
			+ "_s" + std::to_string(m_dictionaryProcessor.sparsity())
			+ "_ox" + std::to_string(m_dictionaryProcessor.patchOffset()[0])
			+ "_p" + std::to_string(m_dictionaryProcessor.patchSize()[0]);
	ASTex::create_directory(ioPath);
	if(m_loadDictionaryMode)
		m_dictionaryProcessor.load(ioPath);
	else
	{
		m_dictionaryProcessor.setInput(m_input);
		m_dictionaryProcessor.dictionaryLearning(m_numberIterationsLearning);
		m_dictionaryProcessor.save(ioPath);
	}
	for(unsigned i=0; i<m_numberResults; ++i)
	{
		ImageType out = m_dictionaryProcessor.synthesize(m_input.width(), m_input.height(), m_numberIterationsSynthesising);
		IO::save01_in_u8(out, ioPath + "/" + "output_it" + std::to_string(m_numberIterationsSynthesising) + "_time" + std::to_string(time(nullptr)) + ".png");
	}
	IO::save01_in_u8(m_input, ioPath + "/" + "input.png");
}


#endif
