.header ON
.mode csv
.output MDD.csv 
-- Extract samples that we will exclude from our analysis
CREATE TEMP TABLE Exclude AS
SELECT DISTINCT sample_id
FROM f20544 mental
WHERE mental.pheno = 2
    OR mental.pheno = 3
    OR mental.pheno = 10;
-- Extract cases
CREATE TEMP TABLE cases AS
SELECT DISTINCT info.sample_id AS sample_id,
    1 AS MDD
FROM (
        SELECT sample_id,
            1 AS pheno,
            0 AS core,
            1 AS essential
        FROM f20436 -- Most of the day or more affected during worst episode of depression
        WHERE pheno > 2
            AND instance = 0 
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id,
            1 AS pheno,
            0 AS core,
            1 AS essential
        FROM f20439 -- Depressed (almost) every day during worst episode of depression
        WHERE pheno > 1
            AND instance = 0 
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id,
            1 AS pheno,
            0 AS core,
            1 AS essential
        FROM f20440 -- More than a little impact on normal roles during worst period of depression
        WHERE pheno > 1
            AND instance = 0 
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id,
            1 AS pheno,
            1 AS core,
            0 AS essential
        FROM f20441 -- Ever had prolonged loss of interest in normal activities
        WHERE pheno = 1
            AND instance = 0 
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id,
            1 AS pheno,
            1 AS core,
            0 AS essential
        FROM f20446 -- Ever had prolonged feelings of sadness or depression
        WHERE pheno = 1
            AND instance = 0 
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id,
            1 AS pheno,
            0 AS core,
            0 AS essential
        FROM f20435 -- Difficulty concentrating during worst episode of depression
        WHERE pheno = 1
            AND instance = 0 
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id,
            1 AS pheno,
            0 AS core,
            0 AS essential
        FROM f20437 -- Thoughts of death during worst episode of depression
        WHERE pheno = 1
            AND instance = 0 
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id,
            1 AS pheno,
            0 AS core,
            0 AS essential
        FROM f20449 -- Feelings of tiredness during worst episode of depression
        WHERE pheno = 1
            AND instance = 0 
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id,
            1 AS pheno,
            0 AS core,
            0 AS essential
        FROM f20450 -- Feelings of worthlessness during worst episode of depression
        WHERE pheno = 1
            AND instance = 0 
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id,
            1 AS pheno,
            0 AS core,
            0 AS essential
        FROM f20532 -- Sleep change during worst episode of depression
        WHERE pheno = 1
            AND instance = 0 
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id,
            1 AS pheno,
            0 AS core,
            0 AS essential
        FROM f20536 -- Weight change during worst episode of depression
        WHERE pheno > 0
            AND instance = 0
    ) AS info
WHERE info.sample_id NOT IN (
        SELECT sample_id
        from Exclude
    )
GROUP BY info.sample_id
HAVING sum(info.pheno) >= 5
    AND sum(info.core) > 0
    AND sum(info.essential) = 3;

-- Symptom of Depression
CREATE TEMP TABLE Symptom AS
SELECT DISTINCT info.sample_id AS sample_id
FROM(
        SELECT sample_id,
            pheno
        FROM f20514 -- Lack of interest or pleasure in doing things
        WHERE pheno > 0 
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id,
            pheno
        FROM f20507 -- Feelings of inadequacy
        WHERE pheno > 0 
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id,
            pheno
        FROM f20510 -- Feelings of depression
        WHERE pheno > 0 
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id,
            pheno
        FROM f20508 -- Trouble concentrating on things
        WHERE pheno > 0 
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id,
            pheno
        FROM f20517 -- Trouble falling or staying asleep, or sleeping too much 
        WHERE pheno > 0 
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id,
            pheno
        FROM f20518 -- Changes in speed or amount of moving or speaking
        WHERE pheno > 0 
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id,
            pheno
        FROM f20519 -- Feelings of tiredness or low energy
        WHERE pheno > 0 
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id,
            pheno
        FROM f20513 -- Thoughts of suicide or self-harm
        WHERE pheno > 0 
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id,
            pheno
        FROM f20511 -- Poor appetite or overeating
        WHERE pheno > 0
    ) AS info
GROUP BY info.sample_id
HAVING sum(info.pheno) > 4;
-- Extract Controls
CREATE TEMP TABLE Controls AS
SELECT s.sample_id AS sample_id,
    0 AS MDD
FROM Participant s
WHERE s.sample_id NOT IN (
        SELECT sample_id
        FROM f20544 -- Report any mental health problems diagnosed by a professional
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id
        FROM cases 
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id
        FROM Exclude 
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id
        FROM Symptom 
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id
        FROM f20002 -- Report depression in previous interview with psychiatric nurse
        WHERE pheno = 1286 -- depression
            OR pheno = 1291 -- mania/bipolar disorder/manic depression
            OR pheno = 1531 -- post-natal depression
            OR pheno = 1289 -- schizophrenia
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id
        FROM f20126 -- Meet previous criteria for depression or bipolar disorder
        WHERE pheno != 0 
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id
        FROM f41202 -- Diagnoses - main ICD10
        WHERE pheno LIKE 'F3__' 
        -------------------------------------------------------------------
        UNION ALL
        SELECT sample_id
        FROM f41204 -- Diagnoses - secondary ICD10
        WHERE pheno LIKE 'F3__'
        UNION ALL
        SELECT sample_id
        FROM f20003 -- Report use of anti-depressant medication at baseline
        WHERE pheno IN (
                1140879616,
                1140921600,
                1140879540,
                1140867878,
                1140916282,
                1140909806,
                1140867888,
                1141152732,
                1141180212,
                1140879634,
                1140867876,
                1140882236,
                1141190158,
                1141200564,
                1140867726,
                1140879620,
                1140867818,
                1140879630,
                1140879628,
                1141151946,
                1140867948,
                1140867624,
                1140867756,
                1140867884,
                1141151978,
                1141152736,
                1141201834,
                1140867690,
                1140867640,
                1140867920,
                1140867850,
                1140879544,
                1141200570,
                1140867934,
                1140867758,
                1140867914,
                1140867820,
                1141151982,
                1140882244,
                1140879556,
                1140867852,
                1140867860,
                1140917460,
                1140867938,
                1140867856,
                1140867922,
                1140910820,
                1140882312,
                1140867944,
                1140867784,
                1140867812,
                1140867668,
                1140867940
            )
    );
INSERT INTO Cases(sample_id, MDD)
SELECT Controls.sample_id AS sample_id,
    Controls.MDD AS MDD
FROM Controls;
CREATE TEMP TABLE Education AS
SELECT edu.sample_id AS sample_id,
    max(
        CASE
            WHEN edu.pheno = 1 THEN 1
            ELSE 0
        END
    ) AS Pheno
FROM f6138 edu
WHERE edu.instance = 0
    AND edu.sample_id NOT IN (
        SELECT sample_id
        FROM f6138
        WHERE instance = 0
            AND pheno = -3
    )
GROUP BY edu.sample_id;
/* Now obtain the required data structure*/
SELECT mdd.sample_id as FID,
    mdd.sample_id as IID,
    sex.Pheno as Sex,
    centre.Pheno as Centre,
    mdd.MDD as MDD,
    bmi.pheno as BMI,
    e.pheno as Education
FROM Cases mdd
    LEFT JOIN f21001 bmi ON mdd.sample_id = bmi.sample_id
    AND bmi.instance = 0
    LEFT JOIN f21003 age ON mdd.sample_id = age.sample_id
    AND age.instance = 0
    LEFT JOIN Education e ON mdd.sample_id = e.sample_id
    JOIN Participant s ON s.sample_id = mdd.sample_id
    AND s.withdrawn = 0
    JOIN f31 sex ON mdd.sample_id = sex.sample_id
    AND sex.instance = 0
    JOIN f54 centre ON mdd.sample_id = centre.sample_id
    AND centre.instance = 0;
.quit
