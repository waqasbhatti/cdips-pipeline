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
       -- add more columns here, including: current LC projectid, field,
       -- datarelease, lcversion, LC path, object tags, variable type, external
       -- tags, external URL, candidate tags, transient tags, comments,
       primary key (objectid)
);


drop table if exists lcinfo;
create table lcinfo (
       version integer not null,
       datarelease integer not null,
       lightcurvecolumns text not null,
       lightcurvesortcol text not null,
       lightcurveapertures text not null,  -- JSON map of aperturename to pixval
       lcbestaperture text,                -- best aperture in
                                           -- lightcurveapertures
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


-- PHOTOMETRY TABLES --


-- using timestamp w/o TZ because of partitioning. see:
-- http://justatheory.com/computers/databases/postgresql/use-timestamptz.html
drop table if exists ism_photometry;
create table ism_photometry (
        -- this MUST be in UTC and will be converted from the JD in FITS headers
       frameutcdt timestamp without time zone not null,
       objectid text not null,
       framekey bigint not null,
       photkey text not null,
       -- image subtraction photometry metadata
       xcc real,                       -- X coordinate on chip
       ycc real,                       -- Y coordinate on chip
       xic real,                       -- X coord on CCD after ISM warp
       yic real,                       -- Y coord on CCD after ISM warp
       bgv real,                       -- background measurement
       bge real,                       -- background measurement error
       fsv real,                       -- source extraction S parameter
       fdv real,                       -- source extraction D parameter
       fkv real,                       -- source extraction K parameter
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
       primary key (frameutcdt, objectid, framekey, photkey)
);


-- using timestamp w/o TZ because of partitioning. see:
-- http://justatheory.com/computers/databases/postgresql/use-timestamptz.html
drop table if exists ap_photometry;
create table ap_photometry (
        -- this MUST be in UTC and will be converted from the JD in FITS headers
       frameutcdt timestamp without time zone not null,
       objectid text not null,
       framekey bigint not null,
       photkey text not null,
       -- this is the photometry data in JSONB format, so we can add columns
       -- later if needed
       photdata jsonb not null,
       primary key (frameutcdt, objectid, framekey, photkey)
);


-- this is the table for calibrated frame information
drop table if exists calibratedframes;
create table calibratedframes (
       framekey bigserial not null,
       entryts timestamp with time zone not null default current_timestamp,
       frameisok bool not null,
       -- frame paths
       fits text not null unique, -- uniqueness enables "on conflict (fits, fistar) do nothing" statements
       fistar text not null unique,
       fiphot text,
       wcs text,
       fitsheader jsonb not null, -- metadata about the image. frame, filter, camera, scope, photometry, environment...
       photinfo jsonb,   -- metadata about preliminary photometry on the image.
       primary key (framekey)
);

create unique index calframes_uindx on calibratedframes (
       framekey, fits
);

-- FIXME: will need to add more columns for aperturephot
drop table if exists photometryinfo;
create table photometryinfo (
       -- this uniquely identifies each type of photometry carried out
       -- this is an MD5 hash of all columns below starting at projectid
       photkey text not null,
       -- metadata
       projectid text not null,
       field text not null,
       photmethod text not null default 'imagesub',
       -- ISM info
       ism_photreftype text not null default 'oneframe',
       ism_photrefframe text not null,
       ism_convsubtype text not null default 'reverse',
       ism_lcapertures text not null,
       ism_convkernel text not null,
       primary key (photkey)
);
