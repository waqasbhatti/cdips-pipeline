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
       projectid integer not null default 0,
       field text not null,
       ccd integer not null,
       isactive bool not null default true,
       unixtime double precision not null,
       framepath text not null,
       jpegpath text not null,
       sval real,
       dval real,
       bgv real,
       ndet integer,
       comment text,
       primary key(field, projectid, ccd, isactive)
);


-- for combined photometric reference selection per field, ccd, and projectid
CREATE TABLE photrefs (
       projectid integer not null default 0,
       field text not null,
       ccd integer not null,
       photreftype text not null,
       isactive integer not null,
       unixtime double precision not null,
       framepath text not null,
       jpegpath text not null,
       convolvetarget text not null,
       convolveregpath text not null,
       cmrawphotpath text not null,
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
