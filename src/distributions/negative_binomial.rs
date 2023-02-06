use num_traits::pow::Pow;
use rand::distributions::Distribution; // for sample
use rand::rngs::StdRng;
use statrs::distribution::Bernoulli;
use statrs::distribution::NegativeBinomial;

pub struct NBinom {
    mu: f64,
    sigma_squared: f64,
    distribution: NegativeBinomial,
}

pub struct ZINB {
    nbinom: NBinom,
    p: f64,
    bernoulli: Bernoulli,
}

fn reparam1(mu: f64, sigma_squared: f64) -> (f64, f64) {
    let mean: f64;
    let variance: f64;

    // safety
    if mu >= sigma_squared {
        mean = mu;
        variance = mu + 1.0;
    } else {
        mean = mu;
        variance = sigma_squared;
    }

    let r: f64 = mean.pow(2) / (variance - mean);
    let p: f64 = mean / variance;

    (r, p)
}

impl NBinom {
    pub fn new(mu: f64, sigma_squared: f64) -> Self {
        let (r, p): (f64, f64) = reparam1(mu, sigma_squared);
        let distribution = NegativeBinomial::new(r, p).unwrap();
        NBinom {
            mu: mu,
            sigma_squared: sigma_squared,
            distribution: distribution,
        }
    }

    pub fn sample(&self, rng: &mut StdRng) -> u64 {
        let count: u64 = self.distribution.sample(rng);
        count
    }
}

impl ZINB {
    pub fn new(mu: f64, sigma_squared: f64, p: f64) -> Self {
        ZINB {
            nbinom: NBinom::new(mu, sigma_squared),
            p: p,
            bernoulli: Bernoulli::new(p).unwrap(),
        }
    }

    pub fn sample(&self, rng: &mut StdRng) -> u64 {
        let count: u64 = if self.bernoulli.sample(rng) == 0.0 {
            self.nbinom.sample(rng)
        } else {
            0
        };

        count
    }
}

impl std::fmt::Display for ZINB {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "(μ: {}, σ²: {}, p: {})",
            self.nbinom.mu, self.nbinom.sigma_squared, self.p
        )
    }
}

impl Clone for NBinom {
    fn clone(&self) -> NBinom {
        NBinom {
            mu: self.mu.clone(),
            sigma_squared: self.sigma_squared.clone(),
            distribution: self.distribution.clone(),
        }
    }
}

impl Clone for ZINB {
    fn clone(&self) -> ZINB {
        ZINB {
            nbinom: self.nbinom.clone(),
            p: self.p.clone(),
            bernoulli: self.bernoulli.clone(),
        }
    }
}
