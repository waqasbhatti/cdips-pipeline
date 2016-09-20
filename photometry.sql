drop table if exists objectinfo;
create table objectinfo (
       objectid text not null,
       twomassid text not null,
       ucac4id text,
       ra double precision not null,
       decl double precision not null,
       pmra double precision,
       pmra_err double precision,
       pmdecl double precision,
       pmdecl_err double precision,
       parallax double precision,
       jmag real,
       hmag real,
       kmag real,
       umag real,
       bmag real,
       vmag real,
       rmag real,
       imag real,
       sdssu real,
       sdssg real,
       sdssr real,
       sdssi real,
       sdssz real,
       gaiag real,
       ndet integer not null,
       network text not null,
       stations text not null,
       primary key (objectid)
);


drop table if exists lcinfo;
create table lcinfo (
       version integer not null,
       datarelease integer not null,
       lightcurvecolumns text not null,
       lightcurvesortcol text not null,
       lightcurveapertures text not null,  -- JSON map of aperturename to pixval
       lcbestaperture text,                -- best aperture in lightcurveapertures
       objectinfocolumns text not null,
       objectinfoidcol text not null,
       unixtime double precision not null,
       lcserver text not null,
       comment text not null,
       primary key (version, datarelease)
);


drop table if exists filters;
create table filters (
       filterid integer not null,
       filtername text not null,
       filterdesc text not null,
       primary key(filterid)
);


drop table if exists lightcurve;
create table lightcurve (
       -- frame metadata
       net text not null,              -- HAT network
       stf integer not null,           -- station number
       prj text not null default '0',  -- project ID of this observation
       fld text not null,              -- observed field ID
       cfn integer not null,           -- camera frame serial number
       cfs text not null default '',   -- camera subframe ID
       ccd integer not null,           -- camera CCD position number
       frt text not null,              -- image frame type (flat, object, etc.)
       -- time
       rjd double precision not null   -- reduced Julian date (midexp)
       bjd double precision,           -- Baryocentric Julian date (midexp)
       -- filter config
       flt integer,                    -- filter ID used in filters table
       flv integer default 0,          -- filter version used
       -- camera config
       cid integer,                    -- camera ID (which camera model?)
       cvn integer default 0,          -- camera version
       cbv integer default 0,          -- camera bias version
       cdv integer default 0,          -- camera dark version
       cfv integer default 0,          -- camera flat version
       exp real,                       -- exposure time in seconds
       -- telescope config
       tid integer,                    -- telescope ID (i.e. which lens?)
       tvn integer default 0,          -- telescope version
       tfs real,                       -- telescope focus setting
       ttt real,                       -- telescope tube temperature
       tms text,                       -- telescope mount state
       tmi integer,                    -- telescope mount ID
       tmv integer default 0,          -- telescope mount version
       tgs text default '',            -- telescope guider status
       -- environment
       mph real,                       -- moonphase at time exposure taken
       iha real,                       -- hour angle of observation
       izd real,                       -- zenith distance of observation
       -- aperture photometry metadata
       xcc real,                       -- X coordinate on chip
       ycc real,                       -- Y coordinate on chip
       bgv real,                       -- background measurement
       bge real,                       -- background measurement error
       fsv real,                       -- source extraction S parameter
       fdv real,                       -- source extraction D parameter
       fkv real,                       -- source extraction K parameter
       -- image subtraction photometry metadata
       ist text not null,              -- image subtype ('reverse'/'normal')
       ipr text not null,              -- ISM photreftype ('one[frame,hour,night]')
       xic real default null,          -- x coord on CCD after ISM warp
       yic real default null           -- y coord on CCD after ISM warp
       -- aperture photometry columns
       aim_000 real,                   -- aperture 000 instrumental mag
       aie_000 real,                   -- aperture 000 instrumental mag err
       aiq_000 text,                   -- aperture 000 instrumental mag flag
       aim_001 real,                   -- aperture 001 instrumental mag
       aie_001 real,                   -- aperture 001 instrumental mag err
       aiq_001 text,                   -- aperture 001 instrumental mag flag
       aim_002 real,                   -- aperture 002 instrumental mag
       aie_002 real,                   -- aperture 002 instrumental mag err
       aiq_002 text,                   -- aperture 002 instrumental mag flag
       arm_000 real,                   -- aperture 000 reduced mag after magfit
       arm_001 real,                   -- aperture 001 reduced mag after magfit
       arm_002 real,                   -- aperture 002 reduced mag after magfit
       aep_000 real,                   -- aperture 000 reduced mag after magfit
       aep_001 real,                   -- aperture 001 reduced mag after magfit
       aep_002 real,                   -- aperture 002 reduced mag after magfit
       atf_000 real,                   -- aperture 000 reduced mag after magfit
       atf_001 real,                   -- aperture 001 reduced mag after magfit
       atf_002 real,                   -- aperture 002 reduced mag after magfit
       -- image subtraction photometry columns
       ifl_000 real,                   -- aperture 000 flux (ADU)
       ife_000 real,                   -- aperture 000 flux err (ADU)
       irm_000 real,                   -- aperture 000 ISM mag (after magfit)
       ire_000 real,                   -- aperture 000 ISM mag err
       irq_000 text,                   -- aperture 000 ISM mag flag
       ifl_001 real,                   -- aperture 001 flux (ADU)
       ife_001 real,                   -- aperture 001 flux err (ADU)
       irm_001 real,                   -- aperture 001 ISM mag (after magfit)
       ire_001 real,                   -- aperture 001 ISM mag err
       irq_001 text,                   -- aperture 001 ISM mag flag
       ifl_002 real,                   -- aperture 002 flux (ADU)
       ife_002 real,                   -- aperture 002 flux err (ADU)
       irm_002 real,                   -- aperture 002 ISM mag (after magfit)
       ire_002 real,                   -- aperture 002 ISM mag err
       irq_002 text,                   -- aperture 002 ISM mag flag
       iep_000 real,                   -- aperture 000 ISM EPD mag
       iep_001 real,                   -- aperture 001 ISM EPD mag
       iep_002 real,                   -- aperture 002 ISM EPD mag
       itf_000 real,                   -- aperture 000 ISM TFA mag
       itf_001 real,                   -- aperture 001 ISM TFA mag
       itf_002 real,                   -- aperture 002 ISM TFA mag
       primary key (net, stf, prj, fld, cfn, cfs, ccd, ist, ipr)
);

-- for each aperture, alter table lightcurve to add the following columns
-- XXX is the aperture index
--       aim_XXX real,                      -- aperture XXX instrumental mag
--       aie_XXX real,                      -- aperture XXX instrumental mag error
--       aiq_XXX text,                      -- aperture XXX instrument mag flag
--       arm_XXX real,                      -- aperture XXX reduced mag after magfit
--       aep_XXX real,                      -- aperture XXX EPD mag
--       atf_XXX real,                      -- aperture XXX TFA mag
--       irm_XXX real,                      -- aperture XXX instrumental mag
--       ire_XXX real,                      -- aperture XXX instrumental mag error
--       irq_XXX text,                      -- aperture XXX instrument mag flag
--       iep_XXX real,                      -- aperture XXX EPD mag
--       itf_XXX real,                      -- aperture XXX TFA mag
