#include "space_weather.hpp"
#include "nrlmsise.hpp"

int dso::NrlmsiseDataHunter::hunt(const dso::MjdEpoch &tt)
{
    /* most probably we are on the same date and 3hour index */
    if ((tt.imjd() == _lastepoch.imjd()) && (_ci == _3hidx(tt)))
    {
        /* data already in _lastdata */
        return 0;
    } /* maybe we are only on the next index ... */
    else if ((tt > _lastepoch) && (tt.diff<DateTimeDifferenceType::FractionalSeconds>(_lastepoch).seconds() < SEC_IN_3H))
    {
        int _3hi = 0;
        if (tt.imjd() != _lastepoch.imjd())
        {
            ++_ci;
        }
        else
        {
            _3hi = _3hidx(tt);
        }
        fill_data(_swdata->cbegin() + _ci, _3hi, &_lastdata);
        /* update last used date */
        _lastepoch = tt;
        return 0;
    } /* or maybe the previous one ...*/
    else if ((tt < _lastepoch) && (_lastepoch.diff<DateTimeDifferenceType::FractionalSeconds>(tt).seconds() < SEC_IN_3H))
    {
        int _3hi = 7;
        if (tt.imjd() != _lastepoch.imjd())
        {
            --_ci;
        }
        else
        {
            _3hi = _3hidx(tt);
        }
        fill_data(_swdata->cbegin() + _ci, _3hi, &_lastdata);
        /* update last used date */
        _lastepoch = tt;
        return 0;
    }
    else
    {
        /* fuck it, search from the top */
        return pr_get_data(tt, &_lastdata);
    }
}

int dso::NrlmsiseDataHunter::fill_data(std::vector<dso::SpaceWeatherData>::const_iterator it, int apidx, dso::Msise00Data *data)
{
    if (!data)
        data = &(this->_lastdata);

    const int _daily_ap_average = dso::SpaceWeatherData::_daily_ap_average;
    const int _3hapindex_start = dso::SpaceWeatherData::_3hapindex_start;

    data->f107A = it->mf107A();
    data->f107 = (it - 1)->mf107();
    data->ap[0] = it->int_array()[_daily_ap_average];
    /* next three are the prior 3hour intervals */
    data->ap[1] = it->int_array()[_3hapindex_start + apidx];
    printf("\t> first index is %d + %d\n", _3hapindex_start, apidx);
    printf("\t> here is the array of ints:");
    for (int i=0; i<22;i++)printf("%d,", it->int_array()[i]);
    printf("\n");
    for (int i = 0; i < 3; i++)
    {
        if (apidx == 0)
        {
            --it;
            apidx = 7;
        }
        else
        {
            --apidx;
        }
        data->ap[1 + i] = it->int_array()[_3hapindex_start + apidx];
    }
    /* next, Average of eight 3 hr AP indicies from 12 to 33 hrs prior to
     * current time */
    double average = 0;
    for (int i = 0; i < 8; i++)
    {
        if (apidx == 0)
        {
            --it;
            apidx = 7;
        }
        else
        {
            --apidx;
        }
        average += it->int_array()[_3hapindex_start + apidx];
    }
    data->ap[5] = average / 8.;
    /* next, Average of eight 3 hr AP indicies from 36 to 57 hrs prior to
     * current time */
    for (int i = 0; i < 8; i++)
    {
        if (apidx == 0)
        {
            --it;
            apidx = 7;
        }
        else
        {
            --apidx;
        }
        average += it->int_array()[_3hapindex_start + apidx];
    }
    data->ap[6] = average / 8.;

    return 0;
}

int dso::NrlmsiseDataHunter::pr_get_data(const dso::MjdEpoch &tt, dso::Msise00Data *mdata)
{
    /* search for day of interest in array */
    for (auto it = _swdata->cbegin(); it != _swdata->cend(); ++it)
    {
        /* if we are on the same MJD, we matched */
        if (tt.imjd() == it->tt().imjd())
        {
            /* we must have data for two days prior to current */
            if (std::distance(_swdata->cbegin(), it) < 2)
            {
                fprintf(stderr,
                        "[ERROR] Not enough space weather data (prior to current "
                        "date) for NRLMSISE00 model! (traceback: %s)\n",
                        __func__);
                return 1;
            }
            /* update current index */
            _ci = std::distance(_swdata->cbegin(), it);
            /* set last epoch used */
            _lastepoch = tt;
            /* fill in values */
            fill_data(it, _3hidx(tt), mdata);
            /* all done */
            return 0;
        }
    }

    /* no date matched! */
    return 1;
}
