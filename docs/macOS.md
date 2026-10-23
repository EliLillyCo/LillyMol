# MacOS

The Mac installation is relatively recent and not completely tested. If you encounter
problems, please raise an issue.

Several software packages are needed.

- install [homebrew](https://brew.sh)
- `brew install bazelisk libomp wget gnutls xz bash gcc`
- `brew install ruby protobuf`
- `gem install google-protobuf`
- clone repo
- `cd /path/to/LillyMol`
- `export LILLYMOL_HOME=/path/to/LillyMol`
- `make -j 10 2>&1 | tee make.log`

bash is needed because the default system bash is too old, but the
version installed by brew is new enough.

Currently the macOS build depends on gcc rather than clang. Work is
underway to address that.

I wish to express my thanks to Harry Stern for providing invaluable help with the
LillyMol MacOS port.

All Bazel unit tests should pass. If you run the tests in the test directory
you may find some "failing" simply due to different orderings of C++
hashes between compiler versions and machines. Eventually we will
address these discrepancies - they are not failures, but poorly
designed tests.

# Resolving Bazel TLS Certificate Errors on Lilly macOS (Zscaler)

## Background

On Lilly-managed Macs, Bazel fails to download dependencies from
`bcr.bazel.build` and `github.com` with errors like:

```
TLS error: (certificate_unknown) PKIX path building failed:
sun.security.provider.certpath.SunCertPathBuilderException:
unable to find valid certification path to requested target
```

This happens because Lilly's network uses **Zscaler**, a corporate TLS
inspection proxy that intercepts HTTPS traffic and re-signs it with a
Zscaler certificate. Bazel runs on an embedded JVM that does not trust
the Zscaler certificate by default.

The fix has three parts:

1. Install a real JDK (Temurin), import the Zscaler certificate into its
   trust store
2. Tell Bazel's embedded JVM to use that trust store for all SSL
   connections
3. Edit .bazelrc to tell bazel how to connect.

---

## Step 1 — Install Eclipse Temurin (OpenJDK)

Bazel bundles its own JVM but has no `keytool` we can use directly.
Install Temurin to get a real `keytool`:

```bash
brew install --cask temurin
```

---

## Step 2 — Extract the Zscaler certificate

Use `echo "Q"` to close the SSL session cleanly so the pipeline
completes:

```bash
echo "Q" | openssl s_client -connect bcr.bazel.build:443 -showcerts 2>/dev/null \
  | awk '/BEGIN CERTIFICATE/,/END CERTIFICATE/{ print }' \
  > /tmp/bcr-chain.pem
```

Inspect the chain to identify the certificates:

```bash
openssl crl2pkcs7 -nocrl -certfile /tmp/bcr-chain.pem \
  | openssl pkcs7 -print_certs -noout
```

Extract the topmost certificate in the chain (the Zscaler intermediate
that acts as the trust anchor in this environment). Replace `N` with the
total number of certificates reported by the `grep` count:

```bash
grep -c "BEGIN CERTIFICATE" /tmp/bcr-chain.pem

# note the count → N
awk '/BEGIN CERTIFICATE/{c++} c==N{ print }' /tmp/bcr-chain.pem \
  > /tmp/zscaler-cert.pem
```

Verify the extracted certificate looks correct:

```bash
openssl x509 -in /tmp/zscaler-cert.pem -noout -subject -issuer
```

---

## Step 3 — Import the certificate into Temurin's trust store

```bash
sudo /Library/Java/JavaVirtualMachines/temurin-26.jdk/Contents/Home/bin/keytool \
  -importcert \
  -alias zscaler-intermediate \
  -file /tmp/zscaler-cert.pem \
  -keystore /Library/Java/JavaVirtualMachines/temurin-26.jdk/Contents/Home/lib/security/cacerts \
  -storepass changeit \
  -noprompt
```

Confirm the import succeeded:

```bash
/Library/Java/JavaVirtualMachines/temurin-26.jdk/Contents/Home/bin/keytool \
  -list \
  -keystore /Library/Java/JavaVirtualMachines/temurin-26.jdk/Contents/Home/lib/security/cacerts \
  -storepass changeit \
  | grep -i "zscaler"
```

You should see a line containing `zscaler-intermediate`.

> **Note:** The Temurin version number (`temurin-26`) may differ if a
> newer version was installed. Adjust the paths above to match the
> version present under `/Library/Java/JavaVirtualMachines/`.

---

## Step 4 — Configure Bazel to use Temurin's trust store

Add the following two lines to your project `.bazelrc` (typically
`LillyMol/src/.bazelrc`):
These lines may already be present, but commented out. Watch for different
versions of temurin.

```
startup --host_jvm_args=-Djavax.net.ssl.trustStore=/Library/Java/JavaVirtualMachines/temurin-26.jdk/Contents/Home/lib/security/cacerts
startup --host_jvm_args=-Djavax.net.ssl.trustStorePassword=changeit
```

This redirects Bazel's embedded JVM to use Temurin's trust store —
where the Zscaler certificate now lives — for all SSL connections.

---

## Step 5 — Restart Bazel and verify

```bash
bazel shutdown
bazel build //...
```

Dependency downloads from `bcr.bazel.build` and `github.com` should now
succeed without TLS errors.
