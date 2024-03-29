# Black-Scholes Model
In order to know more information about a stock option, this options calculator with Black-Scholes Model, the first widely used model for option pricing, can provide the call/put option price, d1, d2, and Greek letters. It can assist investors in establishing an option trading strategy.

## Certain assumptions must be made due to this calculator is modeled by Black-Scholes model.

- It works on European options that can only be exercised at expiration.
- No dividends are paid out during the option’s life.
- Stock markets are efficient. The movement of the markets cannot be predicted, and there is continuous trading
- There are no transaction and commissions costs in buying the option.
- The risk-free rate and volatility of the underlying are known and constant.
- The returns on the underlying are normally distributed.

## Input variables:
- Underlying price (per share)
- Strike price of the option (per share)
- Time to maturity (years)
- Continuously compounding risk-free interest rate (%)
- Volatility (%)

## Output Variables:

The Greek letters, including delta, gamma, vega, rho, and theta, represent the sensitivities of the option price to a single –unit change in the value of either a state variable or a parameter.

The delta of an option is defined as the rate of change of the option price respected to the rate of the change of underlying asset price.
The gamma is defined as the rate of change of delta respected to the rate of change of underlying asset price.
The vega of an option is defined as the rate of change of the option price respected to the volatility of the underlying asset.
The rho of an option is defined as the rate of the option price respected to the interest rate.
The theta of an option is defined as the rate of change of the option price respected to the passage of time.
