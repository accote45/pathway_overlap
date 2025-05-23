.mode csv
.header on 
.output AD.csv

/*
    Get all samples with healthy mum
    If max phenotype is -17, then the mum should not be suffering from 
    other diseases in the list (e.g. heart disease)
    Use the max instance of healthy record, thus provide the latest age of 
    mum
*/
CREATE TEMP TABLE mum_pheno
AS
SELECT mum.sample_id AS sample_id,
            (CASE WHEN (CAST(mum.pheno AS INT))==-17 
            THEN -2
            ELSE 
                CASE WHEN (CAST(mum.pheno AS INT))==10
                THEN 1
                ELSE 
                    CASE WHEN (CAST(mum.pheno AS INT)) BETWEEN 1 AND 9
                        AND ((CAST(mum.pheno AS INT)) BETWEEN 1 AND 9) 
                    THEN 0
                    ELSE -1 /* Don't know or prefer not to answer*/
                    END
                END
            END) AS pheno,
            MIN(mum.instance) AS instance
FROM        f20110 AS mum
WHERE       mum.pheno IN (-17, -11, -13, 1, 2, 6, 8, 9, 10)
GROUP BY mum.sample_id, pheno;

CREATE TEMP TABLE healthy_mum
AS
SELECT  mum.sample_id AS sample_id,
            mum.pheno AS pheno,
            age.pheno AS age,
            death.pheno AS death_age
        FROM (
            SELECT  a.*
            FROM    mum_pheno AS a
            LEFT OUTER JOIN mum_pheno AS b
                ON a.sample_id = b.sample_id AND
                    a.pheno < b.pheno 
            WHERE b.sample_id IS NULL
        ) AS mum
            LEFT JOIN f1845 age ON
                mum.sample_id = age.sample_id
                AND mum.instance = age.instance
            LEFT JOIN f3526 death ON
                mum.sample_id = death.sample_id
                AND mum.instance = death.instance;


/* Repeat the same for father */
CREATE TEMP TABLE dad_pheno
AS
SELECT dad.sample_id AS sample_id,
            (CASE WHEN (CAST(dad.pheno AS INT))==-17 
            THEN -2
            ELSE 
                CASE WHEN (CAST(dad.pheno AS INT))==10
                THEN 1
                ELSE 
                    CASE WHEN (CAST(dad.pheno AS INT)) BETWEEN 1 AND 9
                        AND ((CAST(dad.pheno AS INT)) BETWEEN 1 AND 9) 
                    THEN 0
                    ELSE -1 /* Don't know or prefer not to answer*/
                    END
                END
            END) AS pheno,
            MIN(dad.instance) AS instance
FROM        f20107 AS dad
WHERE       dad.pheno IN (-17, -11, -13, 1, 2, 6, 8, 9, 10)
GROUP BY dad.sample_id, pheno;

CREATE TEMP TABLE healthy_dad
AS
SELECT dad.sample_id AS sample_id,
            dad.pheno AS pheno,
            age.pheno AS age,
            death.pheno AS death_age
        FROM (
            SELECT  a.*
            FROM    dad_pheno AS a
            LEFT OUTER JOIN dad_pheno AS b
                ON a.sample_id = b.sample_id AND
                    a.pheno < b.pheno 
            WHERE b.sample_id IS NULL
        ) AS dad
            LEFT JOIN f2946 age ON
                dad.sample_id = age.sample_id
                AND dad.instance = age.instance
            LEFT JOIN f1807 death ON
                dad.sample_id = death.sample_id
                AND dad.instance = death.instance;
                

/* Exclude samples if they have other neurodegenerative disease that are incompatible with AD */
CREATE  TEMP TABLE Exclude_AD AS
SELECT DISTINCT sample_id 
FROM 
    (
    SELECT  DISTINCT sample_id
    FROM    f41202 
    WHERE   pheno LIKE '"G310"'
            OR pheno LIKE '"F020"'
    UNION 
    SELECT  DISTINCT sample_id   
    FROM    f41204 
    WHERE   pheno LIKE '"G310"'
            OR pheno LIKE '"F020"'
    UNION 
    SELECT  DISTINCT sample_id
    FROM    f40001 
    WHERE   pheno LIKE '"G310"' 
            OR pheno LIKE '"F020"'
    UNION 
    SELECT  DISTINCT sample_id   
    FROM    f40002 
    WHERE   pheno LIKE '"G310"'
            OR pheno LIKE '"F020"'
    UNION 
    SELECT  DISTINCT sample_id
    FROM    gp_clinical 
    WHERE
        read3 IN 
            ('.E115', '.E116',
                '.G78.', 
                'E004.', 'E0040', 'E0041', 'E0042', 'E0043', 'E004z',
                'Eu01.', 'Eu010', 'Eu011', 'Eu012', 'Eu013', 'Eu01y', 'Eu01z',
                'Eu020', 'Eu025',
                'F111.', 'F116.', 'F118.',
                'F11x2', 'F11y2', 'F21y2',
                'X0034', 'X0035', 'X0036', 'X0037', 'X0039', 
                'X003A', 'X003B', 'X003C', 'X003D', 'X003m', 'X003R', 'X003T', 'X003V', 'X003W',
                'Xa0lH', 'Xa0sC', 'Xa0sE', 'XaE74', 'XaKyY', 'XE1Xs')
        OR read2 IN 
            ('E004.', 'E0040', 'E0041', 'E0042', 'E0043', 'E004z',
            'E012.', 'Eu01.', 'Eu010', 'Eu011', 'Eu012', 'Eu013', 'Eu01y', 'Eu01z',
            'Eu020', 'Eu025',
            'F111.', 'F116.', 'F118.',
            'F11x2', 'F11y2', 'F21y2')
    ) as subquery;

/* Obtain the age of diagnosis of all AD cases from GP record*
/* remember to cast the phenotype or we will miss out a lot of entries*/
CREATE      TEMP TABLE AD_GP_Patient AS
SELECT      DISTINCT birth_year.sample_id AS sample_id, 
            MIN(strftime('%Y',gp_clinical.date_event) - CAST(birth_year.pheno AS INT)) AS age,
            1 AS pheno
FROM        f34 birth_year
            LEFT JOIN gp_clinical ON
                gp_clinical.sample_id=birth_year.sample_id
                AND gp_clinical.sample_id IS NOT NULL 
WHERE       (
                gp_clinical.read3 IN 
                    ('.F21Z', 
                        'Eu00.', 'Eu000', 'Eu001', 'Eu002', 'Eu00z', 
                        'F110.', 'F1100', 'F1101', 'Fyu30', 
                        'X002x', 'X002y', 'X002z', 'X0030', 'X0031', 'X0032', 'X0033', 'X003G',
                        'XaIKB', 'XaIKC', 'XE17j')
                OR gp_clinical.read2 IN 
                    ('Eu00.', 'Eu000', 'Eu001', 'Eu002', 'Eu00z', 
                        'F110.', 'F1100', 'F1101', 'Fyu30' )
            )
GROUP BY    birth_year.sample_id
HAVING COUNT(DISTINCT birth_year.pheno)=1; 


/* Get AD patients based on their ICD10 Code*/
CREATE TEMP TABLE AD_ICD_Patient AS
SELECT DISTINCT sample_id,
    1 AS pheno, 
    MIN(age) AS age -- Use the earliest age where AD is diagnosed
FROM (
        SELECT DISTINCT subquery.sample_id AS sample_id, 
                MIN(strftime('%Y', date_record.pheno)-CAST(birth.pheno AS INT)) AS age 
        FROM
        (
            SELECT DISTINCT sample_id, 
                        array 
            FROM        f41202 
            WHERE       pheno LIKE '"G30%' 
                        OR pheno LIKE '"F00%'
            GROUP BY    sample_id
            UNION 
            SELECT DISTINCT sample_id, 
                        array 
            FROM        f41204 
            WHERE       pheno LIKE '"G30%' 
                        OR pheno LIKE '"F00%'
            GROUP BY    sample_id
            UNION
            SELECT DISTINCT sample_id, 
                        array 
            FROM        f41270 
            WHERE       pheno LIKE '"G30%' 
                        OR pheno LIKE '"F00%'
            GROUP BY    sample_id
        ) as subquery   
            LEFT JOIN f41280 date_record ON
                subquery.sample_id = date_record.sample_id 
                AND subquery.array = date_record.array
            LEFT JOIN f34 birth ON 
                subquery.sample_id = birth.sample_id
            GROUP BY subquery.sample_id
    UNION 
    /* Now process death record to obtain samples die of AD and their age */
    SELECT  DISTINCT icd.sample_id AS sample_id,
            MIN(death.pheno) as age
    FROM (
        SELECT DISTINCT sample_id, 
                    instance 
        FROM        f40001 
        WHERE       pheno LIKE '"G30%' 
                    OR pheno LIKE '"F00%'
        GROUP BY    sample_id
        UNION 
        SELECT DISTINCT sample_id, 
                    instance 
        FROM        f40002 
        WHERE       pheno LIKE '"G30%' 
                    OR pheno LIKE '"F00%'
        GROUP BY    sample_id
        ) as icd
        LEFT JOIN f40007 death ON  
            death.sample_id=icd.sample_id 
            AND death.instance=icd.instance
        GROUP BY death.sample_id
        HAVING COUNT( DISTINCT death.pheno)=1
)
GROUP BY sample_id
HAVING COUNT(DISTINCT age)=1;

CREATE  TEMP TABLE AD_Patient AS
SELECT  DISTINCT all_sample.sample_id AS sample_ID,
        all_sample.age as age,
        all_sample.pheno as pheno 
FROM 
(
    SELECT DISTINCT s.sample_id AS sample_id,
                    0 AS pheno, 
                    MAX(CAST(age.pheno AS INT)) AS age /* Want the latest age for individual who haven't developed AD*/
    FROM            Participant AS s
        LEFT JOIN   f21003 age ON
                    s.sample_ID = age.sample_ID
    WHERE s.sample_id NOT IN ( /* Not an AD patient */
        SELECT sample_id FROM AD_ICD_Patient 
        UNION
        SELECT sample_id FROM AD_GP_Patient
    )
    GROUP BY        s.sample_ID
    UNION 
    SELECT  icd.sample_id AS sample_id, 
            icd.pheno AS pheno, 
            icd.age AS age 
            FROM AD_ICD_Patient AS icd
    UNION 
    SELECT  gp.sample_id AS sample_id,
            gp.pheno AS pheno,
            gp.age AS age
            FROM AD_GP_Patient as gp
) AS all_sample
WHERE all_sample.sample_id NOT IN (
    SELECT  sample_id FROM   Exclude_AD 
);

DROP TABLE Exclude_AD;
DROP TABLE AD_ICD_Patient;
DROP TABLE AD_GP_Patient;


SELECT  s.sample_id AS FID,  
        s.sample_id AS IID,
        sex.pheno AS Sex,
        centre.pheno AS Centre,
        mum.pheno AS Mum,
        dad.pheno AS Dad,
        dad.age AS Dad_age,
        mum.age AS Mum_age,
        dad.death_age AS Dad_death,
        mum.death_age AS Mum_death,
        ad.Pheno AS AD,
        SUM(adopted.Pheno) AS Adopted,
        ad.age AS Age
FROM    Participant s 
        LEFT JOIN   f31 sex ON
                    s.sample_id=sex.sample_id 
                    AND sex.instance = 0
        LEFT JOIN   f54 centre ON 
                    s.sample_id=centre.sample_id 
                    AND centre.instance = 0
        LEFT JOIN   healthy_mum mum ON
                    s.sample_id=mum.sample_id
        LEFT JOIN   healthy_dad dad ON
                    s.sample_id=dad.sample_id
        LEFT JOIN   AD_Patient ad ON
                    s.sample_id=ad.sample_id
        LEFT JOIN   f1767 adopted ON 
                    s.sample_id=adopted.sample_id 
                    AND adopted.pheno >=0
group by s.sample_id;
.quit
