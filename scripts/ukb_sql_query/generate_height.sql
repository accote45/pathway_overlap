.mode csv
.header on 
.output Height.csv
SELECT  s.sample_id AS FID, 
        s.sample_id AS IID,
        age.pheno AS Age,
        sex.pheno AS Sex,
        height.pheno AS Height,
        centre.pheno AS Centre,
        ses.pheno AS SES
FROM    Participant s 
        LEFT JOIN   f31 sex ON
                    s.sample_id=sex.sample_id 
                    AND sex.instance = 0
        LEFT JOIN   f21003 age ON 
                    s.sample_id=age.sample_id
                    AND age.instance = 0
        LEFT JOIN   f54 centre ON 
                    s.sample_id=centre.sample_id 
                    AND centre.instance = 0
        LEFT JOIN   f189 ses ON 
                    s.sample_id=ses.sample_id 
                    AND ses.instance = 0
        LEFT JOIN   f50 height ON
                    s.sample_id=height.sample_id 
                    AND height.instance = 0;
.quit
