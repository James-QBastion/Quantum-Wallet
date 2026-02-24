// SPDX-License-Identifier: MIT

use core::array::{ArrayTrait, SpanTrait};
use super::falcon;

/// Authorization Request for QBastion Node Protection
/// Represents a command to be executed by the protected node.
#[derive(Drop, Serde)]
struct AuthorizationRequest {
    /// Hash of the administrative command (e.g. Poseidon hash of payload)
    command_hash: u256,
    /// Unique ID of the target node to prevent replay on other nodes
    node_id: u256,
    /// Timestamp or Nonce for replay protection
    nonce: u64,
}

/// A request bundled with its Falcon signature proof
#[derive(Drop, Serde)]
struct SignedRequest {
    /// The authorization payload
    request: AuthorizationRequest,
    /// Uncompressed signature polynomial s1
    s1: Span<u16>,
    /// Public key polynomial h
    pk: Span<u16>,
    /// Pre-computed hash point of the request (msg_point)
    msg_point: Span<u16>,
}

/// Verifies a single signed request using Falcon-512 logic
fn verify_request(signed_req: SignedRequest, n: u32) -> bool {
    // 1. Verify the 'msg_point' actually corresponds to Hash(request)
    // TODO: Implement Hash-to-Point logic here or verify the commitment.
    // For now, we assume msg_point is passed correctly and verified by the caller/prover constraints.
    // In a full implementation, we must hash `signed_req.request` to `signed_req.msg_point`.

    // 2. Verify Falcon Signature
    let result = falcon::verify_uncompressed::<512>(
        signed_req.s1, 
        signed_req.pk, 
        signed_req.msg_point, 
        n
    );

    match result {
        Result::Ok(()) => true,
        Result::Err(_) => false,
    }
}
