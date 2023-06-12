/*

 This file is part of Butterworth Filter Design, a pair C++ classes and an
 accompanying suite of unit tests for designing high order Butterworth IIR &
 EQ filters using the bilinear transform.
 The generated filter coefficients are split out into cascaded biquad sections,
 for easy use in your garden variety biquad or second-order section (SOS).

 Reference: http://en.wikipedia.org/wiki/Butterworth_filter
 http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt


 Copyright (C) 2013,  iroro orife

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

#include <math.h>

#include "Biquad.h"

Biquad::Biquad()
{
}
Biquad::~Biquad()
{
}

void Biquad::DF2TFourthOrderSection(double B0, double B1, double B2, double B3, double B4,
                                    double A0, double A1, double A2, double A3, double A4)
{
  b0 = B0 / A0;
  b1 = B1 / A0;
  b2 = B2 / A0;
  b3 = B3 / A0;
  b4 = B4 / A0;

  a1 = (-A1) / A0; // The negation conforms the Direct-Form II Transposed discrete-time
  a2 = (-A2) / A0; // filter (DF2T) coefficients to the expectations of the process function.
  a3 = (-A3) / A0;
  a4 = (-A4) / A0;
}

void Biquad::DF2TBiquad(double B0, double B1, double B2,
                        double A0, double A1, double A2)
{
  b0 = B0 / A0;
  b1 = B1 / A0;
  b2 = B2 / A0;
  a1 = (-A1) / A0; // The negation conforms the Direct-Form II Transposed discrete-time
  a2 = (-A2) / A0; // filter (DF2T) coefficients to the expectations of the process function.
}

void BiquadChain::allocate(uint32_t count)
{

  numFilters = count;
  m_yn.resize(numFilters);
  m_yn1.resize(numFilters);
  m_yn2.resize(numFilters);

  // Fourth order sections
  m_yn3.resize(numFilters);
  m_yn4.resize(numFilters);
}

BiquadChain::BiquadChain() : numFilters(0)
{
}

BiquadChain::BiquadChain(uint32_t count)
{
  allocate(count);
  reset();
}
BiquadChain::~BiquadChain()
{
}

void BiquadChain::resize(uint32_t count)
{
  allocate(count);
}

void BiquadChain::reset()
{

  m_xn1 = 0;
  m_xn2 = 0;

  for (uint32_t i = 0; i < numFilters; i++)
  {
    m_yn[i] = 0;
    m_yn1[i] = 0;
    m_yn2[i] = 0;

    // Fourth order sections
    m_yn3[i] = 0;
    m_yn4[i] = 0;
  }
}

void BiquadChain::processBiquad(const double *input, double *output, const int32_t stride, const uint32_t count, const Biquad *coeffs)
{

  double *yn = &m_yn[0];
  double *yn1 = &m_yn1[0];
  double *yn2 = &m_yn2[0];

  for (uint32_t n = 0; n < count; n++)
  {
    double xn = *input;

    yn[0] = coeffs[0].b0 * xn +
            coeffs[0].b1 * m_xn1 +
            coeffs[0].b2 * m_xn2 +
            coeffs[0].a1 * yn1[0] +
            coeffs[0].a2 * yn2[0];

    for (uint32_t i = 1; i < numFilters; i++)
    {
      yn[i] = coeffs[i].b0 * yn[i - 1] +
              coeffs[i].b1 * yn1[i - 1] +
              coeffs[i].b2 * yn2[i - 1] +
              coeffs[i].a1 * yn1[i] +
              coeffs[i].a2 * yn2[i];
    }

    // Shift delay line elements.
    for (uint32_t i = 0; i < numFilters; i++)
    {
      yn2[i] = yn1[i];
      yn1[i] = yn[i];
    }
    m_xn2 = m_xn1;
    m_xn1 = xn;

    // Store result and stride
    *output = yn[numFilters - 1];

    input += stride;
    output += stride;
  }
}

void BiquadChain::processFourthOrderSections(const double *input, double *output, const int32_t stride, const uint32_t count, const Biquad *coeffs)
{

  double *yn = &m_yn[0];
  double *yn1 = &m_yn1[0];
  double *yn2 = &m_yn2[0];
  double *yn3 = &m_yn3[0];
  double *yn4 = &m_yn4[0];

  for (uint32_t n = 0; n < count; n++)
  {
    double xn = *input;

    yn[0] = coeffs[0].b0 * xn +
            coeffs[0].b1 * m_xn1 +
            coeffs[0].b2 * m_xn2 +
            coeffs[0].b3 * m_xn3 +
            coeffs[0].b4 * m_xn4 +

            coeffs[0].a1 * yn1[0] +
            coeffs[0].a2 * yn2[0] +
            coeffs[0].a3 * yn3[0] +
            coeffs[0].a4 * yn4[0];

    for (uint32_t i = 1; i < numFilters; i++)
    {
      yn[i] = coeffs[i].b0 * yn[i - 1] +
              coeffs[i].b1 * yn1[i - 1] +
              coeffs[i].b2 * yn2[i - 1] +
              coeffs[i].b3 * yn3[i - 1] +
              coeffs[i].b4 * yn4[i - 1] +

              coeffs[i].a1 * yn1[i] +
              coeffs[i].a2 * yn2[i] +
              coeffs[i].a3 * yn3[i] +
              coeffs[i].a4 * yn4[i];
    }

    // Shift delay line elements.
    for (uint32_t i = 0; i < numFilters; i++)
    {
      yn4[i] = yn3[i];
      yn3[i] = yn2[i];
      yn2[i] = yn1[i];
      yn1[i] = yn[i];
    }

    m_xn4 = m_xn3;
    m_xn3 = m_xn2;
    m_xn2 = m_xn1;
    m_xn1 = xn;

    // Store result and stride
    *output = yn[numFilters - 1];

    input += stride;
    output += stride;
  }
}
