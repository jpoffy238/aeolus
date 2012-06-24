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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "rankwave.h"

extern float exp2ap(float);

Rngen Pipewave::_rgen;
float *Pipewave::_arg = 0;
float *Pipewave::_att = 0;

void Pipewave::initstatic(float fsamp) {
	int k;

	if (_arg)
		return;
	k = (int) (fsamp);
	_arg = new float[k];
	k = (int) (0.5f * fsamp);
	_att = new float[k];
}

void Pipewave::play(void) {
	int i, d, k1, k2;
	float g, dg, y, dy, t;
	float *p, *q, *r;

	p = play_pointer;
	r = release_pointer;

	if (_sdel & 1) {
		if (!p) {
			p = attack_start;
			play_interpolation = 0.0f;
			play_interpolation_speed = 0.0f;
		}
	} else {
		if (!r) {
			r = p;
			p = 0;
			release_gain = 1.0f;
			release_interpolation = play_interpolation;
			release_count = release_length;
		}
	}

	if (r) {
		k1 = PERIOD;
		q = _out;
		g = release_gain;
		i = release_count - 1;
		dg = g / PERIOD;
		if (i)
			dg *= release_multiplier;

		if (r < loop_start) {
			while (k1--) {
				*q++ += g * *r++;
				g -= dg;
			}
		} else {
			y = release_interpolation;
			dy = release_detune;
			while (k1) {
				t = y + k1 * dy;
				d = 0;
				k2 = k1;
				if (t > 1.0f) {
					d = 1;
					k2 = (int) ((1.0f - y) / dy);
				} else if (t < 0.0f) {
					d = -1;
					k2 = (int) (-y / dy);
				}
				k1 -= k2;
				while (k2--) {
					*q++ += g * (r[0] + y * (r[1] - r[0]));
					g -= dg;
					y += dy;
					r += sample_step;
				}
				y -= d;
				r += d;
			}
			release_interpolation = y;
		}

		if (i) {
			release_gain = g;
			release_count = i;
			if (r >= loop_end)
				r -= loop_length;
		} else
			r = 0;
	}

	if (p) {
		k1 = PERIOD;
		q = _out;
		if (p < loop_start) {
			while (k1--) {
				*q++ += *p++;
			}
		} else {
			y = play_interpolation;
			play_interpolation_speed += instability * 0.0005f
					* (0.05f * instability * (_rgen.urandf() - 0.5f) - play_interpolation_speed);
			dy = play_interpolation_speed * sample_step;
			while (k1) {
				t = y + k1 * dy;
				d = 0;
				k2 = k1;
				if (t > 1.0f) {
					d = 1;
					k2 = (int) ((1.0f - y) / dy);
				} else if (t < 0.0f) {
					d = -1;
					k2 = (int) (-y / dy);
				}
				k1 -= k2;
				while (k2--) {
					*q++ += p[0] + y * (p[1] - p[0]);
					y += dy;
					p += sample_step;
				}
				y -= d;
				p += d;
			}
			if (p >= loop_end)
				p -= loop_length;
			play_interpolation = y;
		}
	}

	play_pointer = p;
	release_pointer = r;
}

void Pipewave::generateWaves(AdditiveSynth *D, int n, float fsamp, float fpipe) {
	int h, i, k, nc;
	float f0, f1, f, m, t, v, v0;

	m = D->_n_att.vi(n);
	for (h = 0; h < NUMBER_OF_HARMONICS; h++) {
		t = D->_h_att.vi(h, n);
		if (t > m)
			m = t;
	}
	attack_length = (int) (fsamp * m + 0.5);
	attack_length = (attack_length + PERIOD - 1) & ~(PERIOD - 1);

	f1 = (fpipe + D->_n_off.vi(n) + D->_n_ran.vi(n) * (2 * _rgen.urand() - 1))
			/ fsamp;
	f0 = f1 * exp2ap(D->_n_atd.vi(n) / 1200.0f);

	for (h = NUMBER_OF_HARMONICS - 1; h >= 0; h--) {
		f = (h + 1) * f1;
		if ((f < 0.45f) && (D->_h_lev.vi(h, n) >= -40.0f))
			break;
	}
	if (f > 0.250f)
		sample_step = 3;
	else if (f > 0.125f)
		sample_step = 2;
	else
		sample_step = 1;

	looplen(f1 * fsamp, sample_step * fsamp, (int) (fsamp / 6.0f), &loop_length, &nc);
	if (loop_length < sample_step * PERIOD) {
		k = (sample_step * PERIOD - 1) / loop_length + 1;
		loop_length *= k;
		nc *= k;
	}

	k = attack_length + loop_length + sample_step * (PERIOD + 4);

	delete[] attack_start;
	attack_start = new float[k];
	loop_start = attack_start + attack_length;
	loop_end = loop_start + loop_length;
	memset(attack_start, 0, k * sizeof(float));

	release_length = (int) (ceil(D->_n_dct.vi(n) * fsamp / PERIOD) + 1);
	release_multiplier = 1.0f - pow(0.1, 1.0 / release_length);
	release_detune = sample_step * (exp2ap(D->_n_dcd.vi(n) / 1200.0f) - 1.0f);
	instability = D->_n_ins.vi(n);

	t = 0.0f;
	k = (int) (fsamp * D->_n_att.vi(n) + 0.5);
	for (i = 0; i <= attack_length; i++) {
		_arg[i] = t - floor(t + 0.5);
		t += (i < k) ? (((k - i) * f0 + i * f1) / k) : f1;
	}

	for (i = 1; i < loop_length; i++) {
		t = _arg[attack_length] + (float) i * nc / loop_length;
		_arg[i + attack_length] = t - floor(t + 0.5);
	}

	v0 = exp2ap(0.1661 * D->_n_vol.vi(n));
	for (h = 0; h < NUMBER_OF_HARMONICS; h++) {
		if ((h + 1) * f1 > 0.45)
			break;
		v = D->_h_lev.vi(h, n);
		if (v < -80.0)
			continue;

		v = v0
				* exp2ap(
						0.1661
								* (v
										+ D->_h_ran.vi(h, n)
												* (2 * _rgen.urand() - 1)));
		k = (int) (fsamp * D->_h_att.vi(h, n) + 0.5);
		attgain(k, D->_h_atp.vi(h, n));

		for (i = 0; i < attack_length + loop_length; i++) {
			t = _arg[i] * (h + 1);
			t -= floor(t);
			m = v * sin(2 * M_PI * t);
			if (i < k)
				m *= _att[i];
			attack_start[i] += m;
		}
	}
	for (i = 0; i < sample_step * (PERIOD + 4); i++)
		attack_start[i + attack_length + loop_length] = attack_start[i + attack_length];
}

void Pipewave::looplen(float f, float fsamp, int lmax, int *aa, int *bb) {
	int i, j, a, b, t;
	int z[8];
	double g, d;

	g = fsamp / f;
	for (i = 0; i < 8; i++) {
		a = z[i] = (int) (floor(g + 0.5));
		g -= a;
		b = 1;
		j = i;
		while (j > 0) {
			t = a;
			a = z[--j] * a + b;
			b = t;
		}
		if (a < 0) {
			a = -a;
			b = -b;
		}
		if (a <= lmax) {
			d = fsamp * b / a - f;
			if ((fabs(d) < 0.1) && (fabs(d) < 3e-4 * f))
				break;
			g = (fabs(g) < 1e-6) ? 1e6 : 1.0 / g;
		} else {
			b = (int) (lmax * f / fsamp);
			a = (int) (b * fsamp / f + 0.5);
			d = fsamp * b / a - f;
			break;
		}
	}
	*aa = a;
	*bb = b;
}

void Pipewave::attgain(int n, float p) {
	int i, j, k;
	float d, m, w, x, y, z;

	w = 0.05;
	x = 0.0;
	y = 0.6;
	if (p > 0)
		y += 0.11 * p;
	z = 0.0;
	j = 0;
	for (i = 1; i <= 24; i++) {
		k = n * i / 24;
		x = 1.0 - z - 1.5 * y;
		y += w * x;
		d = w * y * p / (k - j);
		while (j < k) {
			m = (double) j / n;
			_att[j++] = (1.0 - m) * z + m;
			z += d;
		}
	}
}

void Pipewave::save(FILE *F) {
	int k;
	union {
		int16_t i16[16];
		int32_t i32[8];
		float flt[8];
	} d;

	d.i32[0] = attack_length;
	d.i32[1] = loop_length;
	d.i16[4] = sample_step;
	d.i16[5] = release_length;
	d.flt[3] = release_multiplier;
	d.i32[4] = 0;
	d.i32[5] = 0;
	d.i32[6] = 0;
	d.i32[7] = 0;
	fwrite(&d, 1, 32, F);
	k = attack_length + loop_length + sample_step * (PERIOD + 4);
	fwrite(attack_start, k, sizeof(float), F);
}

void Pipewave::load(FILE *F) {
	int k;
	union {
		int16_t i16[16];
		int32_t i32[8];
		float flt[8];
	} d;

	fread(&d, 1, 32, F);
	attack_length = d.i32[0];
	loop_length = d.i32[1];
	sample_step = d.i16[4];
	release_length = d.i16[5];
	release_multiplier = d.flt[3];
	k = attack_length + loop_length + sample_step * (PERIOD + 4);
	delete[] attack_start;
	attack_start = new float[k];
	loop_start = attack_start + attack_length;
	loop_end = loop_start + loop_length;
	fread(attack_start, k, sizeof(float), F);
}

Rankwave::Rankwave(int n0, int n1) :
		_n0(n0), _n1(n1), _list(0), _modif(false) {
	_pipes = new Pipewave[n1 - n0 + 1];
}

Rankwave::~Rankwave(void) {
	delete[] _pipes;
}

void Rankwave::generateWaves(AdditiveSynth *additiveSynth, float fsamp, float fbase, float *scale) {
	Pipewave::initstatic(fsamp);

	fbase *= additiveSynth->_fn / (additiveSynth->_fd * scale[9]);
	for (int i = _n0; i <= _n1; i++) {
		_pipes[i - _n0].generateWaves(additiveSynth, i - _n0, fsamp, ldexp(fbase * scale[i % 12], i / 12 - 5));
	}
	_modif = true;
}

void Rankwave::set_param(float *out, int del, int pan) {
	int n, a, b;
	Pipewave *P;

	_sbit = 1 << del;
	switch (pan) {
	case 'L':
		a = 2, b = 0;
		break;
	case 'C':
		a = 2, b = 1;
		break;
	case 'R':
		a = 2, b = 2;
		break;
	default:
		a = 4, b = 0;
	}
	for (n = _n0, P = _pipes; n <= _n1; n++, P++)
		P->_out = out + ((n % a) + b) * PERIOD;
}

void Rankwave::play(int shift) {
	Pipewave *P, *Q;

	for (P = 0, Q = _list; Q; Q = Q->_link) {
		Q->play();
		if (shift)
			Q->_sdel = (Q->_sdel >> 1) | Q->_sbit;
		if (Q->_sdel || Q->play_pointer || Q->release_pointer)
			P = Q;
		else {
			if (P)
				P->_link = Q->_link;
			else
				_list = Q->_link;
		}
	}
}

int Rankwave::save(const char *path, AdditiveSynth *D, float fsamp, float fbase,
		float *scale) {
	FILE *F;
	Pipewave *P;
	int i;
	char name[1024];
	char data[64];
	char *p;

	sprintf(name, "%s/%s", path, D->_filename);
	if ((p = strrchr(name, '.')))
		strcpy(p, ".ae1");
	else
		strcat(name, ".ae1");

	F = fopen(name, "wb");
	if (F == NULL) {
		fprintf(stderr, "Can't open waveform file '%s' for writing\n", name);
		return 1;
	}

	memset(data, 0, 16);
	strcpy(data, "ae1");
	data[4] = 1;
	fwrite(data, 1, 16, F);

	memset(data, 0, 64);
	data[0] = 0;
	data[1] = 0;
	data[2] = 0;
	data[3] = 0;
	data[4] = _n0;
	data[5] = _n1;
	data[6] = 0;
	data[7] = 0;
	*((float *) (data + 8)) = fsamp;
	*((float *) (data + 12)) = fbase;
	memcpy(data + 16, scale, 12 * sizeof(float));
	fwrite(data, 1, 64, F);

	for (i = _n0, P = _pipes; i <= _n1; i++, P++)
		P->save(F);

	fclose(F);

	_modif = false;
	return 0;
}

int Rankwave::load(const char *path, AdditiveSynth *D, float fsamp, float fbase,
		float *scale) {
	FILE *F;
	Pipewave *P;
	int i;
	char name[1024];
	char data[64];
	char *p;
	float f;

	sprintf(name, "%s/%s", path, D->_filename);
	if ((p = strrchr(name, '.')))
		strcpy(p, ".ae1");
	else
		strcat(name, ".ae1");

	F = fopen(name, "rb");
	if (F == NULL) {
#ifdef DEBUG
		fprintf (stderr, "Can't open waveform file '%s' for reading\n", name);
#endif
		return 1;
	}

	fread(data, 1, 16, F);
	if (strcmp(data, "ae1")) {
#ifdef DEBUG
		fprintf (stderr, "File '%s' is not an Aeolus waveform file\n", name);
#endif
		fclose(F);
		return 1;
	}

	if (data[4] != 1) {
#ifdef DEBUG
		fprintf (stderr, "File '%s' has an incompatible version tag (%d)\n", name, data [4]);
#endif
		fclose(F);
		return 1;
	}

	fread(data, 1, 64, F);
	if (_n0 != data[4] || _n1 != data[5]) {
#ifdef DEBUG
		fprintf (stderr, "File '%s' has an incompatible note range (%d %d), (%d %d)\n", name, _n0, _n1, data [4], data [5]);
#endif
		fclose(F);
		return 1;
	}

	f = *((float *) (data + 8));
	if (fabsf(f - fsamp) > 0.1f) {
#ifdef DEBUG
		fprintf (stderr, "File '%s' has a different sample frequency (%3.1lf)\n", name, f);
#endif
		fclose(F);
		return 1;
	}

	f = *((float *) (data + 12));
	if (fabsf(f - fbase) > 0.1f) {
#ifdef DEBUG
		fprintf (stderr, "File '%s' has a different tuning (%3.1lf)\n", name, f);
#endif
		fclose(F);
		return 1;
	}

	for (i = 0; i < 12; i++) {
		f = *((float *) (data + 16 + 4 * i));
		if (fabsf(f / scale[i] - 1.0f) > 6e-5f) {
#ifdef DEBUG
			fprintf (stderr, "File '%s' has a different temperament\n", name);
#endif
			fclose(F);
			return 1;
		}
	}

	for (i = _n0, P = _pipes; i <= _n1; i++, P++)
		P->load(F);

	fclose(F);

	_modif = false;
	return 0;
}
