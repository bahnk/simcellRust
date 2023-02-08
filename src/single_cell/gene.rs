use crate::distributions::negative_binomial::ZINB;

pub struct Gene {
    pub id: String,
    pub symbol: String,
    pub status: String,
    pub distribution: ZINB,
}

impl Gene {
    pub fn new(id: String, symbol: String, status: String, distribution: ZINB) -> Self {
        Gene {
            id: id,
            symbol: symbol,
            status: status,
            distribution: distribution,
        }
    }
}

impl std::fmt::Display for Gene {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "({}, {}, {}, {})",
            self.id, self.symbol, self.status, self.distribution
        )
    }
}

impl Clone for Gene {
    fn clone(&self) -> Gene {
        Gene {
            id: self.id.clone(),
            symbol: self.symbol.clone(),
            status: self.status.clone(),
            distribution: self.distribution.clone(),
        }
    }
}
