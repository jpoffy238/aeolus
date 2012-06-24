/*
 Copyright (C) 2003-2010 Fons Adriaensen <fons@kokkinizita.net>

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

#include <sys/time.h>
#include "clalsadrv.h"

// Public members ----------------------------------------------------------------------

Alsa_driver::~Alsa_driver(void) {
	snd_pcm_sw_params_free(_capt_swpar);
	snd_pcm_hw_params_free(_capt_hwpar);
	snd_pcm_sw_params_free(_play_swpar);
	snd_pcm_hw_params_free(_play_hwpar);

	if (_play_handle)
		snd_pcm_close(_play_handle);
	if (_capt_handle)
		snd_pcm_close(_capt_handle);
	if (_ctrl_handle)
		snd_ctl_close(_ctrl_handle);
}

Alsa_driver::Alsa_driver(const char *play_name, const char *capt_name,
		const char *ctrl_name, unsigned int rate, snd_pcm_uframes_t frsize,
		unsigned int nfrags) :
		_rate(rate), _frsize(frsize), _nfrags(nfrags), _play_handle(0), _capt_handle(
				0), _ctrl_handle(0), _play_hwpar(0), _play_swpar(0), _capt_hwpar(
				0), _capt_swpar(0), _play_nchan(0), _capt_nchan(0), _play_npfd(
				0), _capt_npfd(0), _synced(false), _stat(-1) {
	initialise(play_name, capt_name, ctrl_name, rate, frsize, nfrags);
}

Alsa_driver::Alsa_driver(const char *name, unsigned int rate,
		snd_pcm_uframes_t frsize, unsigned int nfrags, bool play, bool capt,
		bool ctrl) :
		_rate(rate), _frsize(frsize), _nfrags(nfrags), _play_handle(0), _capt_handle(
				0), _ctrl_handle(0), _play_hwpar(0), _play_swpar(0), _capt_hwpar(
				0), _capt_swpar(0), _play_nchan(0), _capt_nchan(0), _play_npfd(
				0), _capt_npfd(0), _synced(false), _stat(-1) {
	initialise(play ? name : 0, capt ? name : 0, ctrl ? name : 0, rate, frsize,
			nfrags);
}

int Alsa_driver::pcm_start(void) {
	unsigned int i, j, n;
	int err;

	if (_play_handle) {
		n = snd_pcm_avail_update(_play_handle);
		if (n != _frsize * _nfrags) {
			fprintf(stderr,
					"Alsa_driver: full buffer not available at start.\n");
			return -1;
		}

		for (i = 0; i < _nfrags; i++) {
			play_init(_frsize);
			for (j = 0; j < _play_nchan; j++)
				clear_chan(j, _frsize);
			play_done(_frsize);
		}

		if ((err = snd_pcm_start(_play_handle)) < 0) {
			fprintf(stderr, "Alsa_driver: pcm_start(play): %s.\n",
					snd_strerror(err));
			return -1;
		}
	}

	if (_capt_handle && !_synced && ((err = snd_pcm_start(_capt_handle)) < 0)) {
		fprintf(stderr, "Alsa_driver: pcm_start(capt): %s.\n",
				snd_strerror(err));
		return -1;
	}

	return 0;
}

int Alsa_driver::pcm_stop(void) {
	int err;

	if (_play_handle && ((err = snd_pcm_drop(_play_handle)) < 0)) {
		fprintf(stderr, "Alsa_driver: pcm_drop(play): %s.\n",
				snd_strerror(err));
		return -1;
	}

	if (_capt_handle && !_synced && ((err = snd_pcm_drop(_capt_handle)) < 0)) {
		fprintf(stderr, "Alsa_driver: pcm_drop(capt): %s.\n",
				snd_strerror(err));
		return -1;
	}

	return 0;
}

snd_pcm_sframes_t Alsa_driver::pcm_wait(void) {
	bool need_capt;
	bool need_play;
	snd_pcm_sframes_t capt_av;
	snd_pcm_sframes_t play_av;
	unsigned short rev;
	int i, r, n1, n2;

	_stat = 0;
	need_capt = _capt_handle ? true : false;
	need_play = _play_handle ? true : false;

	while (need_play || need_capt) {
		n1 = 0;
		if (need_play) {
			snd_pcm_poll_descriptors(_play_handle, _pfd, _play_npfd);
			n1 += _play_npfd;
		}
		n2 = n1;
		if (need_capt) {
			snd_pcm_poll_descriptors(_capt_handle, _pfd + n1, _capt_npfd);
			n2 += _capt_npfd;
		}

		for (i = 0; i < n2; i++)
			_pfd[i].events |= POLLERR;

		r = poll(_pfd, n2, 1000);
		if (r < 0) {
			if (errno == EINTR) {
				_stat = -1;
				return 0;
			}
			fprintf(stderr, "Alsa_driver: poll(): %s\n.", strerror(errno));
			_stat = -2;
			return 0;
		}
		if (r == 0) {
			fprintf(stderr, "Alsa_driver: poll timed out\n.");
			_stat = -3;
			return 0;
		}

		if (need_play) {
			snd_pcm_poll_descriptors_revents(_play_handle, _pfd, n1, &rev);
			if (rev & POLLERR) {
				fprintf(stderr, "Alsa_driver: error on playback pollfd.\n");
				_stat = -4;
				recover();
				return 0;
			}
			if (rev & POLLOUT)
				need_play = false;
		}
		if (need_capt) {
			snd_pcm_poll_descriptors_revents(_capt_handle, _pfd + n1, n2 - n1,
					&rev);
			if (rev & POLLERR) {
				fprintf(stderr, "Alsa_driver: error on capture pollfd.\n");
				_stat = -4;
				recover();
				return 0;
			}
			if (rev & POLLIN)
				need_capt = false;
		}
	}

	play_av = 999999;
	if (_play_handle && (play_av = snd_pcm_avail_update(_play_handle)) < 0) {
		_stat = -5;
		recover();
		return 0;
	}

	capt_av = 999999;
	if (_capt_handle && (capt_av = snd_pcm_avail_update(_capt_handle)) < 0) {
		_stat = -6;
		recover();
		return 0;
	}

	return (capt_av < play_av) ? capt_av : play_av;
}

int Alsa_driver::pcm_idle(snd_pcm_uframes_t len) {
	unsigned int i;
	snd_pcm_uframes_t n, k;

	if (_capt_handle) {
		n = len;
		while (n) {
			k = capt_init(n);
			capt_done(k);
			n -= k;
		}
	}

	if (_play_handle) {
		n = len;
		while (n) {
			k = play_init(n);
			for (i = 0; i < _play_nchan; i++)
				clear_chan(i, k);
			play_done(k);
			n -= k;
		}
	}

	return 0;
}

int Alsa_driver::play_init(snd_pcm_uframes_t len) {
	unsigned int i;
	const snd_pcm_channel_area_t *a;
	int err;

	if ((err = snd_pcm_mmap_begin(_play_handle, &a, &_play_offs, &len)) < 0) {
		fprintf(stderr, "Alsa_driver: snd_pcm_mmap_begin(play): %s.\n",
				snd_strerror(err));
		return -1;
	}

	_play_step = (a->step) >> 3;
	for (i = 0; i < _play_nchan; i++, a++) {
		_play_ptr[i] = (char *) a->addr
				+ ((a->first + a->step * _play_offs) >> 3);
	}

	return len;
}

int Alsa_driver::capt_init(snd_pcm_uframes_t len) {
	unsigned int i;
	const snd_pcm_channel_area_t *a;
	int err;

	if ((err = snd_pcm_mmap_begin(_capt_handle, &a, &_capt_offs, &len)) < 0) {
		fprintf(stderr, "Alsa_driver: snd_pcm_mmap_begin(capt): %s.\n",
				snd_strerror(err));
		return -1;
	}

	_capt_step = (a->step) >> 3;
	for (i = 0; i < _capt_nchan; i++, a++) {
		_capt_ptr[i] = (char *) a->addr
				+ ((a->first + a->step * _capt_offs) >> 3);
	}

	return len;
}

void Alsa_driver::printinfo(void) {
	fprintf(stderr, "playback :");
	if (_play_handle) {
		fprintf(stderr, "\n  nchan  : %d\n", _play_nchan);
		fprintf(stderr, "  rate   : %d\n", _rate);
		fprintf(stderr, "  frsize : %ld\n", _frsize);
		fprintf(stderr, "  nfrags : %d\n", _nfrags);
		fprintf(stderr, "  format : %s\n", snd_pcm_format_name(_play_format));
	} else
		fprintf(stderr, " not enabled\n");
	fprintf(stderr, "capture  :");
	if (_capt_handle) {
		fprintf(stderr, "\n  nchan  : %d\n", _capt_nchan);
		fprintf(stderr, "  rate   : %d\n", _rate);
		fprintf(stderr, "  frsize : %ld\n", _frsize);
		fprintf(stderr, "  nfrags : %d\n", _nfrags);
		fprintf(stderr, "  format : %s\n", snd_pcm_format_name(_capt_format));
		if (_play_handle)
			fprintf(stderr, "%s\n", _synced ? "synced" : "not synced");
	} else
		fprintf(stderr, " not enabled\n");
}

// Private members ---------------------------------------------------------------------

void Alsa_driver::initialise(const char *play_name, const char *capt_name,
		const char *ctrl_name, unsigned int rate, snd_pcm_uframes_t frsize,
		unsigned int nfrags) {
	int err, dir;
	snd_ctl_card_info_t *card;

	if (play_name) {
		if (snd_pcm_open(&_play_handle, play_name, SND_PCM_STREAM_PLAYBACK, 0)
				< 0) {
			_play_handle = 0;
			fprintf(stderr,
					"Alsa_driver: Cannot open PCM device %s for playback.\n",
					play_name);
		}
	}

	if (capt_name) {
		if (snd_pcm_open(&_capt_handle, capt_name, SND_PCM_STREAM_CAPTURE, 0)
				< 0) {
			_capt_handle = 0;
			fprintf(stderr,
					"Alsa_driver: Cannot open PCM device %s for capture.\n",
					capt_name);
		}
	}

	if (!_play_handle && !_capt_handle)
		return;

	if (ctrl_name) {
		snd_ctl_card_info_alloca(&card);

		if ((err = snd_ctl_open(&_ctrl_handle, ctrl_name, 0)) < 0) {
			fprintf(stderr, "Alse_driver: ctl_open(): %s\n", snd_strerror(err));
			return;
		}

		if ((err = snd_ctl_card_info(_ctrl_handle, card)) < 0) {
			fprintf(stderr, "Alsa_driver: ctl_card_info(): %s\n",
					snd_strerror(err));
			return;
		}
	}

	// check capabilities here

	if (_play_handle) {
		if (snd_pcm_hw_params_malloc(&_play_hwpar) < 0) {
			fprintf(stderr, "Alsa_driver: can't allocate playback hw params\n");
			return;
		}

		if (snd_pcm_sw_params_malloc(&_play_swpar) < 0) {
			fprintf(stderr, "Alsa_driver: can't allocate playback sw params\n");
			return;
		}

		if (set_hwpar(_play_handle, _play_hwpar, "playback", &_play_nchan) < 0)
			return;

		if (_play_nchan > MAXPLAY) {
			fprintf(stderr,
					"Alsa_driver: detected %d playback channels, reset to %d.\n",
					_play_nchan, MAXPLAY);
			_play_nchan = MAXPLAY;
		}

		if (set_swpar(_play_handle, _play_swpar, "playback") < 0)
			return;
	}

	if (_capt_handle) {
		if (snd_pcm_hw_params_malloc(&_capt_hwpar) < 0) {
			fprintf(stderr, "Alsa_driver: can't allocate capture hw params\n");
			return;
		}

		if (snd_pcm_sw_params_malloc(&_capt_swpar) < 0) {
			fprintf(stderr, "Alsa_driver: can't allocate capture sw params\n");
			return;
		}

		if (set_hwpar(_capt_handle, _capt_hwpar, "capture", &_capt_nchan) < 0)
			return;

		if (_capt_nchan > MAXCAPT) {
			fprintf(stderr,
					"Alsa_driver: detected %d capture channels, reset to %d\n",
					_capt_nchan, MAXCAPT);
			_capt_nchan = MAXCAPT;
		}

		if (set_swpar(_capt_handle, _capt_swpar, "capture") < 0)
			return;
	}

	if (_play_handle) {
		if (snd_pcm_hw_params_get_rate(_play_hwpar, &rate, &dir)
				|| (rate != _rate) || dir) {
			fprintf(stderr,
					"Alsa_driver: can't get requested sample rate for playback.\n");
			return;
		}

		if (snd_pcm_hw_params_get_period_size(_play_hwpar, &frsize, &dir)
				|| (frsize != _frsize) || dir) {
			fprintf(stderr,
					"Alsa_driver: can't get requested period size for playback.\n");
			return;
		}

		if (snd_pcm_hw_params_get_periods(_play_hwpar, &nfrags, &dir)
				|| (nfrags != _nfrags) || dir) {
			fprintf(stderr,
					"Alsa_driver: can't get requested number of periods for playback.\n");
			return;
		}

		snd_pcm_hw_params_get_format(_play_hwpar, &_play_format);
		snd_pcm_hw_params_get_access(_play_hwpar, &_play_access);

		switch (_play_format) {
		case SND_PCM_FORMAT_S32_LE:
			_clear_func = clear_32le;
			_play_func = play_32le;
			break;

		case SND_PCM_FORMAT_S24_3LE:
			_clear_func = clear_24le;
			_play_func = play_24le;
			break;

		case SND_PCM_FORMAT_S16_LE:
			_clear_func = clear_16le;
			_play_func = play_16le;
			break;

		default:
			fprintf(stderr,
					"Alsa_driver: can't handle playback sample format.\n");
			return;
		}

		_play_npfd = snd_pcm_poll_descriptors_count(_play_handle);
	}

	if (_capt_handle) {
		if (snd_pcm_hw_params_get_rate(_capt_hwpar, &rate, &dir)
				|| (rate != _rate) || dir) {
			fprintf(stderr,
					"Alsa_driver: can't get requested sample rate for capture.\n");
			return;
		}

		if (snd_pcm_hw_params_get_period_size(_capt_hwpar, &frsize, &dir)
				|| (frsize != _frsize) || dir) {
			fprintf(stderr,
					"Alsa_driver: can't get requested period size for capture.\n");
			return;
		}

		if (snd_pcm_hw_params_get_periods(_capt_hwpar, &nfrags, &dir)
				|| (nfrags != _nfrags) || dir) {
			fprintf(stderr,
					"Alsa_driver: can't get requested number of periods for capture.\n");
			return;
		}

		if (_play_handle)
			_synced = !snd_pcm_link(_play_handle, _capt_handle);

		snd_pcm_hw_params_get_format(_capt_hwpar, &_capt_format);
		snd_pcm_hw_params_get_access(_capt_hwpar, &_capt_access);

		switch (_capt_format) {
		case SND_PCM_FORMAT_S32_LE:
			_capt_func = capt_32le;
			break;

		case SND_PCM_FORMAT_S24_3LE:
			_capt_func = capt_24le;
			break;

		case SND_PCM_FORMAT_S16_LE:
			_capt_func = capt_16le;
			break;

		default:
			fprintf(stderr,
					"Alsa_driver: can't handle capture sample format.\n");
			return;
		}

		_capt_npfd = snd_pcm_poll_descriptors_count(_capt_handle);
	}

	if (_play_npfd + _capt_npfd > MAXPFD) {
		fprintf(stderr, "Alsa_driver: interface requires more than %d pollfd\n",
				MAXPFD);
		return;
	}

	_stat = 0;
}

int Alsa_driver::set_hwpar(snd_pcm_t *handle, snd_pcm_hw_params_t *hwpar,
		const char *sname, unsigned int *nchan) {
	int err;
	unsigned int n;

	if ((err = snd_pcm_hw_params_any(handle, hwpar)) < 0) {
		fprintf(stderr, "Alsa_driver: no %s hw configurations available: %s.\n",
				sname, snd_strerror(err));
		return -1;
	}

	if ((err = snd_pcm_hw_params_set_periods_integer(handle, hwpar)) < 0) {
		fprintf(stderr,
				"Alsa_driver: can't set %s period size to integral value.\n",
				sname);
		return -1;
	}

	if (((err = snd_pcm_hw_params_set_access(handle, hwpar,
			SND_PCM_ACCESS_MMAP_NONINTERLEAVED)) < 0)
			&& ((err = snd_pcm_hw_params_set_access(handle, hwpar,
					SND_PCM_ACCESS_MMAP_INTERLEAVED)) < 0)
			&& ((err = snd_pcm_hw_params_set_access(handle, hwpar,
					SND_PCM_ACCESS_MMAP_COMPLEX)) < 0)) {
		fprintf(stderr,
				"Alsa_driver: the %s interface doesn't support mmap-based access.\n",
				sname);
		return -1;
	}

	if (((err = snd_pcm_hw_params_set_format(handle, hwpar,
			SND_PCM_FORMAT_S32_LE)) < 0)
			&& ((err = snd_pcm_hw_params_set_format(handle, hwpar,
					SND_PCM_FORMAT_S24_3LE)) < 0)
			&& ((err = snd_pcm_hw_params_set_format(handle, hwpar,
					SND_PCM_FORMAT_S16_LE)) < 0)) {
		fprintf(stderr,
				"Alsa_driver: the %s interface doesn't support 32, 24 or 16 bit access.\n.",
				sname);
		return -1;
	}

	if ((err = snd_pcm_hw_params_set_rate(handle, hwpar, _rate, 0)) < 0) {
		fprintf(stderr, "Alsa_driver: can't set %s sample rate to %u.\n", sname,
				_rate);
		return -1;
	}

	snd_pcm_hw_params_get_channels_max(hwpar, nchan);

	if (*nchan > 1024) {
		fprintf(stderr,
				"Alsa_driver: detected more than 1024 %s channnels, reset to 2.\n",
				sname);
		*nchan = 2;
	}

	if ((err = snd_pcm_hw_params_set_channels(handle, hwpar, *nchan)) < 0) {
		fprintf(stderr, "Alsa_driver: can't set %s channel count to %u.\n",
				sname, *nchan);
		return -1;
	}

	if ((err = snd_pcm_hw_params_set_period_size(handle, hwpar, _frsize, 0))
			< 0) {
		fprintf(stderr, "Alsa_driver: can't set %s period size to %lu.\n",
				sname, _frsize);
		return -1;
	}

	// This is to handle recent ALSA releases creating a default device with
	// a large number of periods...
	n = _nfrags;
	snd_pcm_hw_params_set_periods_min(handle, hwpar, &n, NULL);
	if (_nfrags < n)
		_nfrags = n;

	if ((err = snd_pcm_hw_params_set_periods_near(handle, hwpar, &_nfrags, 0))
			< 0) {
		fprintf(stderr, "Alsa_driver: can't set %s periods to %u.\n", sname,
				_nfrags);
		return -1;
	}

	if ((err = snd_pcm_hw_params_set_buffer_size(handle, hwpar,
			_frsize * _nfrags)) < 0) {
		fprintf(stderr, "Alsa_driver: can't set %s buffer length to %lu.\n",
				sname, _frsize * _nfrags);
		return -1;
	}

	if ((err = snd_pcm_hw_params(handle, hwpar)) < 0) {
		fprintf(stderr, "Alsa_driver: can't set %s hardware parameters.\n",
				sname);
		return -1;
	}

	return 0;
}

int Alsa_driver::set_swpar(snd_pcm_t *handle, snd_pcm_sw_params_t *swpar,
		const char *sname) {
	int err;

	snd_pcm_sw_params_current(handle, swpar);

	if ((err = snd_pcm_sw_params_set_tstamp_mode(handle, swpar,
			SND_PCM_TSTAMP_MMAP)) < 0) {
		fprintf(stderr, "Alsa_driver: can't set %s timestamp mode to %u.\n",
				sname, SND_PCM_TSTAMP_MMAP);
		return -1;
	}

	if ((err = snd_pcm_sw_params_set_avail_min(handle, swpar, _frsize)) < 0) {
		fprintf(stderr, "Alsa_driver: can't set %s availmin to %lu.\n", sname,
				_frsize);
		return -1;
	}

	if ((err = snd_pcm_sw_params(handle, swpar)) < 0) {
		fprintf(stderr, "Alsa_driver: can't set %s software parameters.\n",
				sname);
		return -1;
	}

	return 0;
}

int Alsa_driver::recover(void) {
	int err;
	snd_pcm_status_t *stat;

	snd_pcm_status_alloca(&stat);

	if ((err = snd_pcm_status(_play_handle ? _play_handle : _capt_handle, stat))
			< 0) {
		fprintf(stderr, "Alsa_driver: pcm_status(): %s\n", snd_strerror(err));
	} else if (snd_pcm_status_get_state(stat) == SND_PCM_STATE_XRUN) {
		struct timeval tnow, trig;

		gettimeofday(&tnow, 0);
		snd_pcm_status_get_trigger_tstamp(stat, &trig);
		fprintf(stderr,
				"Alsa_driver: stat = %02x, xrun of at least %8.3lf ms\n", _stat,
				1e3 * tnow.tv_sec - 1e3 * trig.tv_sec + 1e-3 * tnow.tv_usec
						- 1e-3 * trig.tv_usec);
	}

	if (pcm_stop())
		return -1;

	if (_play_handle && ((err = snd_pcm_prepare(_play_handle)) < 0)) {
		fprintf(stderr, "Alsa_driver: pcm_prepare(play): %s\n",
				snd_strerror(err));
		return -1;
	}

	if (_capt_handle && !_synced
			&& ((err = snd_pcm_prepare(_capt_handle)) < 0)) {
		fprintf(stderr, "Alsa_driver: pcm_prepare(capt): %s\n",
				snd_strerror(err));
		return -1;
	}

	if (pcm_start())
		return -1;

	return 0;
}

// Static members ----------------------------------------------------------------------

char *Alsa_driver::clear_16le(char *dst, int step, int nfrm) {
	while (nfrm--) {
		*((short int *) dst) = 0;
		dst += step;
	}
	return dst;
}

char *Alsa_driver::play_16le(const float *src, char *dst, int step, int nfrm) {
	float s;
	short int d;

	while (nfrm--) {
		s = *src++;
		if (s > 1)
			d = 0x7fff;
		else if (s < -1)
			d = 0x8001;
		else
			d = (short int) (0x7fff * s);
		*((short int *) dst) = d;
		dst += step;
	}
	return dst;
}

const char *Alsa_driver::capt_16le(const char *src, float *dst, int step,
		int nfrm) {
	float d;
	short int s;

	while (nfrm--) {
		s = *((short int *) src);
		d = (float) s / 0x7fff;
		*dst++ = d;
		src += step;
	}
	return src;
}

char *Alsa_driver::clear_24le(char *dst, int step, int nfrm) {
	while (nfrm--) {
		dst[0] = 0;
		dst[1] = 0;
		dst[2] = 0;
		dst += step;
	}
	return dst;
}

char *Alsa_driver::play_24le(const float *src, char *dst, int step, int nfrm) {
	float s;
	int d;

	while (nfrm--) {
		s = *src++;
		if (s > 1)
			d = 0x007fffff;
		else if (s < -1)
			d = 0x00800001;
		else
			d = (int) (0x007fffff * s);
		dst[0] = d;
		dst[1] = d >> 8;
		dst[2] = d >> 16;
		dst += step;
	}
	return dst;
}

const char *Alsa_driver::capt_24le(const char *src, float *dst, int step,
		int nfrm) {
	float d;
	int s;

	while (nfrm--) {
		s = src[0] & 0xFF;
		s += (src[1] & 0xFF) << 8;
		s += (src[2] & 0xFF) << 16;
		if (s & 0x00800000)
			s -= 0x01000000;
		d = (float) s / 0x007fffff;
		*dst++ = d;
		src += step;
	}
	return src;
}

char *Alsa_driver::clear_32le(char *dst, int step, int nfrm) {
	while (nfrm--) {
		*((int *) dst) = 0;
		dst += step;
	}
	return dst;
}

char *Alsa_driver::play_32le(const float *src, char *dst, int step, int nfrm) {
	float s;
	int d;

	while (nfrm--) {
		s = *src++;
		if (s > 1)
			d = 0x007fffff;
		else if (s < -1)
			d = 0x00800001;
		else
			d = (int) (0x007fffff * s);
		*((int *) dst) = d << 8;
		dst += step;
	}
	return dst;
}

const char *Alsa_driver::capt_32le(const char *src, float *dst, int step,
		int nfrm) {
	float d;
	int s;

	while (nfrm--) {
		s = *((int *) src);
		d = (float) s / 0x7fffff00;
		*dst++ = d;
		src += step;
	}
	return src;
}

