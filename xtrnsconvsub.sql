-- this is the table for arefshifted frame information
drop table if exists arefshiftedframes;
create table arefshiftedframes (
       framekey bigserial not null,
       entryts timestamp with time zone not null default current_timestamp,
       -- original fits, astromref, itrans used, output xtrns fits
       arefshiftedframe text not null,
       origframekey bigint not null,
       astomref text not null,
       itrans text not null,
       shiftisok bool not null,
       -- arefshifted frame info
       warpcheckmargin real,
       warpcheckthresh real,
       warpinfopickle bytea,
       primary key (framekey)
);

create unique index arefshiftedframes_uindx on arefshiftedframes (
       framekey, origframekey, astromref, itrans, arefshiftedframe
);


-- this is the table for subtracted frame information
drop table if exists subtractedframes;
create table subtractedframes (
       framekey bigserial not null,
       entryts timestamp with time zone not null default current_timestamp,
       -- original fits, photref, output kernel, output subtraced fits
       subtractedframe text not null,
       origframekey bigint not null,
       arefshiftedkey bigint not null,
       photreftype text not null,
       photref text not null,
       kernel text not null,
       subisok bool not null,
       -- subtracted frame info
       convkernelspec text not null,
       primary key (framekey)
);

create unique index subtractedframes_uindx on subtractedframes (
       framekey, origframekey, astromref, itrans, subtractedframe
);
