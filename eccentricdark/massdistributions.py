# The functions here define the binary mass distribution functions...

import numpy as np

def mass_distribution_sampler(
    form,
    args=None):

    if form=='flat':
        m1 = np.random.uniform(
            low=args[0],
            high=args[1],
            size=1
        ) 
        m2 = np.random.uniform(
            low=args[0],
            high=args[1],
            size=1
        )
        return (m1, m2) #Units M_solar

    elif form=='gaussian':
        m1 = np.random.normal(
            loc=args[0],
            scale=args[1],
            size=1
        ) 
        m2 = np.random.normal(
            loc=args[0],
            scale=args[1],
            size=1
        )
        return (m1, m2) #Units M_solar      

    elif form=='1602.03842': #Double check this 
        def m1_pdf(m1): #To do: not used... 
            if m1<5:
                return 0.
            elif m1>95:
                return 0.
            elif ((m1>=5) and (m1<=95)):
                normalization = (
                    1.35
                    / (np.power(5., -1.35) - np.power(95., -1.35))
                )
                return (normalization * np.power(m1, -2.35))
            else:
                return False
        def m1_cdf(m1): #To do: not used...
            if m1<5:
                return 0.
            elif m1>95:
                return 1.
            elif ((m1>=5) and (m1<=95)):
                return (
                    (np.power(5., -1.35) - np.power(m1, -1.35))
                    / (np.power(5., -1.35) - np.power(95., -1.35))
                )
            else:
                return False
        def m1_invcdf(cdf):
            return (
                np.power(
                    np.power(5., -1.35)
                    - cdf * (
                        np.power(5.,-1.35)
                        - np.power(95., -1.35)
                    ),
                -1./1.35)
            )

        random_sample = np.random.random(1)
        m1 = m1_invcdf(random_sample)

        qmin = (5./m1)
        qmax = 1.

        def q_pdf(q): #To do: not used... 
            if (q<(5./m1)):
                return 0.
            elif (q>1.):
                return 0
            elif ((q>=(5./m1)) and (q<=1.)):
                return 1./(qmax-qmin)
            else:
                return False

        def q_cdf(q): #To do: not used...
            if (q<(5./m1)):
                return 0.
            elif (q>1.):
                return 1
            elif ((q>=(5./m1)) and (q<=1.)):
                return (1./(qmax-qmin))*(q-qmin)
            else:
                return False

        def q_invcdf(cdf):
            return (cdf*(qmax-qmin) + qmin)

        random_sample = np.random.random(1)
        q = q_invcdf(random_sample)
        m2 = m1 * q

        return (m1, m2)

    elif form=='1907.02283':
        def m1_pdf(m1):
            if m1<5:
                return 0.
            elif m1>50:
                return 0.
            elif ((m1>=5) and (m1<=50)):
                normalization = (
                    1.3 /
                    (np.power(5., -1.3) - np.power(50., -1.3))
                )
                return (normalization * np.power(m1, -2.3))
            else:
                return False

        def m1_cdf(m1):
            if m1<5:
                return 0.
            elif m1>50:
                return 1.
            elif ((m1>=5) and (m1<=50)):
                return (
                    (np.power(5., -1.3) - np.power(m1, -1.3))
                    / (np.power(5., -1.3) - np.power(50., -1.3))
                )
            else:
                return False

        def m1_invcdf(cdf):
            return np.power(
                np.power(5., -1.3)
                - cdf * (
                    np.power(5.,-1.3)
                    - np.power(50., -1.3)
                ), -1./1.3
            )

        random_sample = np.random.random(1)
        m1  = m1_invcdf(random_sample)
        m2 = float(m1)

        return (m1, m2)
                    
