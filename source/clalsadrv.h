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

#ifndef __CLALSADRV_H
#define __CLALSADRV_H

#define ALSA_PCM_NEW_HW_PARAMS_API
#define ALSA_PCM_NEW_SW_PARAMS_API

#include <alsa/asoundlib.h>

class Alsa_driver {
public:

	Alsa_driver(const char *play_name, const char *capt_name,
			const char *ctrl_name, unsigned int rate, snd_pcm_uframes_t frsize,
			unsigned int nfrags);

	Alsa_driver(const char *name, unsigned int rate, snd_pcm_uframes_t frsize,
			unsigned int nfrags, bool play, bool capt, bool ctrl);

	~Alsa_driver(void);

	void printinfo(void);

	int pcm_start(void);

	int pcm_stop(void);

	snd_pcm_sframes_t pcm_wait(void);

	int pcm_idle(snd_pcm_uframes_t len);

	int play_init(snd_pcm_uframes_t len);

	void play_chan(int chan, const float *src, snd_pcm_uframes_t len) {
		_play_ptr[chan] = _play_func(src, _play_ptr[chan], _play_step, len);
	}

	void clear_chan(int chan, snd_pcm_uframes_t len) {
		_play_ptr[chan] = _clear_func(_play_ptr[chan], _play_step, len);
	}

	int play_done(snd_pcm_uframes_t len) {
		return snd_pcm_mmap_commit(_play_handle, _play_offs, len);
	}

	int capt_init(snd_pcm_uframes_t len);

	void capt_chan(int chan, float *dst, snd_pcm_uframes_t len) {
		_capt_ptr[chan] = _capt_func(_capt_ptr[chan], dst, _capt_step, len);
	}

	int capt_done(snd_pcm_uframes_t len) {
		return snd_pcm_mmap_commit(_capt_handle, _capt_offs, len);
	}

	int stat(void) {
		return _stat;
	}
	int nplay(void) {
		return _play_nchan;
	}
	int ncapt(void) {
		return _capt_nchan;
	}
	snd_pcm_t *play_handle(void) {
		return _play_handle;
	}
	snd_pcm_t *capt_handle(void) {
		return _capt_handle;
	}

private:

	typedef char *(*clear_function)(char *, int, int);
	typedef char *(*play_function)(const float *, char *, int, int);
	typedef const char *(*capt_function)(const char *, float *, int, int);

	enum {
		MAXPFD = 16, MAXPLAY = 64, MAXCAPT = 64
	};

	void initialise(const char *play_name, const char *capt_name,
			const char *ctrl_name, unsigned int rate, snd_pcm_uframes_t frsize,
			unsigned int nfrags);

	int set_hwpar(snd_pcm_t *handle, snd_pcm_hw_params_t *hwpar,
			const char *sname, unsigned int *nchan);
	int set_swpar(snd_pcm_t *handle, snd_pcm_sw_params_t *swpar,
			const char *sname);
	int recover(void);

	unsigned int _rate;
	snd_pcm_uframes_t _frsize;
	unsigned int _nfrags;
	snd_pcm_t *_play_handle;
	snd_pcm_t *_capt_handle;
	snd_ctl_t *_ctrl_handle;
	snd_pcm_hw_params_t *_play_hwpar;
	snd_pcm_sw_params_t *_play_swpar;
	snd_pcm_hw_params_t *_capt_hwpar;
	snd_pcm_sw_params_t *_capt_swpar;
	snd_pcm_format_t _play_format;
	snd_pcm_format_t _capt_format;
	snd_pcm_access_t _play_access;
	snd_pcm_access_t _capt_access;
	unsigned int _play_nchan;
	unsigned int _capt_nchan;
	int _play_npfd;
	int _capt_npfd;
	bool _synced;
	struct pollfd _pfd[MAXPFD];
	snd_pcm_uframes_t _capt_offs;
	snd_pcm_uframes_t _play_offs;
	int _play_step;
	int _capt_step;
	char *_play_ptr[MAXPLAY];
	const char *_capt_ptr[MAXCAPT];
	int _stat;
	int _pcnt;
	bool _xrun;
	clear_function _clear_func;
	play_function _play_func;
	capt_function _capt_func;

	static char *clear_32le(char *dst, int step, int nfrm);
	static char *clear_24le(char *dst, int step, int nfrm);
	static char *clear_16le(char *dst, int step, int nfrm);
	static char *play_32le(const float *src, char *dst, int step, int nfrm);
	static char *play_24le(const float *src, char *dst, int step, int nfrm);
	static char *play_16le(const float *src, char *dst, int step, int nfrm);
	static const char *capt_32le(const char *src, float *dst, int step,
			int nfrm);
	static const char *capt_24le(const char *src, float *dst, int step,
			int nfrm);
	static const char *capt_16le(const char *src, float *dst, int step,
			int nfrm);
};

#endif

