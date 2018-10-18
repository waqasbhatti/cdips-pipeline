-- PURPOSE: make the TM-refinfo.sqlite DB used by autoimagesub.get_astromref
-- WB notes (2018/10/17) that though postgres can be not, it's not worth
-- depricating this quite yet!

-- database to use for automatic ISM photometry

-- for astrometric reference selection per field, ccd, and projectid
CREATE TABLE astromrefs (
       field text,
       projectid integer,
       ccd integer,
       isactive integer,
       unixtime double precision,
       framepath text,
       jpegpath text,
       sval real,
       dval real,
       bgv real,
       ndet integer,
       comment text,
       primary key(field, projectid, ccd, isactive)
);


-- for combined photometric reference selection per field, ccd, and projectid
CREATE TABLE photrefs (
       field text,
       projectid integer,
       ccd integer,
       photreftype text,
       isactive integer,
       unixtime double precision,
       framepath text,
       jpegpath text,
       convolvetarget text,
       convolveregpath text,
       cmrawphotpath text,
       target_zenithdist real,
       target_moondist real,
       target_moonelev real,
       target_moonphase real,
       target_hourangle real,
       target_ndet real,
       target_medmagerr real,
       target_magerrmad real,
       target_medsrcbgv real,
       target_stdsrcbgv real,
       target_medsval real,
       target_meddval real,
       photrefinfo text, -- photref selection-info pickle JSON
       primary key(field, projectid, ccd, photreftype, isactive)
);
