-- this contains the PostgreSQL DB photometry tables for image subtraction

-- the iphotfile table
create table iphotfiles (
       projectid integer not null default 0,
       field text not null,
       ccd integer not null,
       isactive bool not null default true,
       iphotfilepath text not null,
       framerjd double precision not null,
       framefilepath text not null,
       entrytimestamp timestamp with time zone not null default current_timestamp,
       primary key (projectid, field, ccd, isactive, iphotfilepath)
);


-- the objects table
create table iphotobjects (
       projectid integer not null,
       field text not null,
       ccd integer not null,
       isactive bool not null default true,
       objectid text not null,
       iphotfilepath text not null,
       iphotfileline integer not null,
       entrytimestamp timestamp with time zone not null default current_timestamp,
       primary key (projectid, field, ccd, isactive,
                    objectid, iphotfilepath, iphotfileline)
);


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
       photrefinfo jsonb, -- photref selection-info pickle JSON
       primary key(field, projectid, ccd, photreftype, isactive)
);
