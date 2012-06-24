/*
 Copyright (C) 2003-2008 Fons Adriaensen <fons@kokkinizita.net>

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef __RANKWAVE_H
#define __RANKWAVE_H

#include "addsynth.h"
#include "rngen.h"

#define PERIOD 64

class Pipewave {
private:

	Pipewave(void) :
			attack_start(0), loop_start(0), loop_end(0), loop_length(0), sample_step(0), release_length(0), release_multiplier(0), _link(0), _sbit(
					0), _sdel(0), play_pointer(0), play_interpolation(0), play_interpolation_speed(0), release_pointer(0), release_interpolation(0), release_gain(
					0), release_count(0) {
	}

	~Pipewave(void) {
		delete[] attack_start;
	}

	friend class Rankwave;

	void generateWaves(AdditiveSynth *D, int n, float fsamp, float fpipe);
	void save(FILE *F);
	void load(FILE *F);
	void play(void);

	static void looplen(float f, float fsamp, int lmax, int *aa, int *bb);
	static void attgain(int n, float p);

	float *attack_start; // attack start
	float *loop_start; // loop start
	float *loop_end; // loop end
	int32_t attack_length; // attack length
	int32_t loop_length; // loop length
	int16_t sample_step; // sample step
	int16_t release_length; // release lenght
	float release_multiplier; // release multiplier
	float release_detune; // release detune
	float instability; // instability

	Pipewave *_link; // link to next in active chain
	uint32_t _sbit; // on state bit
	uint32_t _sdel; // delayed state
	float *_out; // audio output buffer
	float *play_pointer; // play pointer
	float play_interpolation; // play interpolation
	float play_interpolation_speed; // play interpolation speed
	float *release_pointer; // release pointer
	float release_interpolation; // release interpolation
	float release_gain; // release gain
	int16_t release_count; // release count

	static void initstatic(float fsamp);

	static Rngen _rgen;
	static float *_arg;
	static float *_att;
};

class Rankwave {
public:

	Rankwave(int n0, int n1);
	~Rankwave(void);

	void note_on(int n) {
		if ((n < _n0) || (n > _n1))
			return;
		Pipewave *P = _pipes + (n - _n0);
		P->_sbit = _sbit;
		if (!(P->_sdel || P->play_pointer || P->release_pointer)) {
			P->_sdel |= _sbit;
			P->_link = _list;
			_list = P;
		}
	}

	void note_off(int n) {
		if ((n < _n0) || (n > _n1))
			return;
		Pipewave *P = _pipes + (n - _n0);
		P->_sdel >>= 4;
		P->_sbit = 0;
	}

	void all_off(void) {
		Pipewave *P;
		for (P = _list; P; P = P->_link)
			P->_sbit = 0;
	}

	int n0(void) const {
		return _n0;
	}
	int n1(void) const {
		return _n1;
	}
	void play(int shift);
	void set_param(float *out, int del, int pan);
	void generateWaves(AdditiveSynth *D, float fsamp, float fbase, float *scale);
	int save(const char *path, AdditiveSynth *D, float fsamp, float fbase,
			float *scale);
	int load(const char *path, AdditiveSynth *D, float fsamp, float fbase,
			float *scale);
	bool modif(void) const {
		return _modif;
	}

	int _cmask; // used by division logic
	int _nmask; // used by division logic

private:

	Rankwave(const Rankwave&);
	Rankwave& operator=(const Rankwave&);

	int _n0;
	int _n1;
	uint32_t _sbit;
	Pipewave *_list;
	Pipewave *_pipes;
	bool _modif;
};

#endif

