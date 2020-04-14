#ifndef __COLOR_MATCHING_H__
#define __COLOR_MATCHING_H__

#include "histogram.h"

template<typename I>
class ColorMatching
{
public:
	ColorMatching();

	//mandatory

	void setInput(const I& input);
	void setModel(const I& model);

	I generate();

private:

	I m_input;
	I m_model;
};

template<typename I>
void ColorMatching<I>::setInput(const I& input)
{
	m_input = input;
}

template<typename I>
void ColorMatching<I>::setModel(const I& model)
{
	m_model = model;
}

template<typename I>
I ColorMatching<I>::generate()
{
}

#endif
